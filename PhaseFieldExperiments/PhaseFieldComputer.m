classdef PhaseFieldComputer < handle

    properties (Constant, Access = public)
        tolErrU = 1e-8;
        tolErrPhi = 1e-8;
        tolErrStag = 1e-8;
        tau = 100*1e2;
    end

    properties (Access = private)
        mesh
        boundaryConditions
        initPhi
        functional
    end

    properties (Access = private)
        shallPrint
        monitor
        data
    end

    methods (Access = public)

        function obj = PhaseFieldComputer(cParams)
            obj.init(cParams)
        end

        function output = compute(obj)
            u = LagrangianFunction.create(obj.mesh,2,'P1');
            phi = obj.initPhi;
            uOld = u;
            phiOld = phi;
          
            iter = 0;
            costFun = null(2,1);
            output = [];

            maxSteps = length(obj.boundaryConditions.bcValues);
            for i = 1:maxSteps
                fprintf('\n ********* STEP %i *********  \n',i)
                bc = obj.boundaryConditions.nextStep();
                u.fValues = obj.updateInitialDisplacement(bc,uOld);

                eStag = 1; iterStag = 1; costOldStag = 0;
                iterUMax = 1; iterPhiMax = 1;
                while (eStag > obj.tolErrStag) && (iterStag < 300)

                    eU = 1; iterU = 1; costOldU = 0;
                    while (eU > obj.tolErrU) && (iterU < 100)
                        LHS = obj.functional.computeElasticLHS(u,phi);
                        RHS = obj.functional.computeElasticRHS(u,phi,bc);
                        u.fValues = obj.computeDisplacement(LHS,RHS,u,bc);
                        F = obj.computeForceVector(LHS,u);

                        [eU, costU] = obj.computeErrorCostFunction(u,phi,bc,costOldU);
                        costOldU = costU;
                        obj.printCost('iterU',iterU,costU,eU);
                        iterU = iterU + 1;

                        costFun(1,end+1) = costU;
                        costFun(2,end) = 0;
                        iter = iter+1;
                        obj.monitor.update(iter,{[],[],[],[],[costFun(1,end)],[]});

                    end
                    if iterU > iterUMax
                        iterUMax = iterU;
                    end
                   
                    ePhi = 1;  iterPhi = 1; costOldPhi = 0;
                    while (ePhi > obj.tolErrPhi) && (iterPhi < 100)
                        LHS = obj.functional.computePhaseFieldLHS(u,phi);
                        RHS = obj.functional.computePhaseFieldRHS(u,phi);
                        phiNew = obj.updatePhi(LHS, RHS, phi.fValues);
                        phi.fValues = obj.projectInLowerAndUpperBound(phiNew,phiOld.fValues,1);

                        [ePhi, costPhi] = obj.computeErrorCostFunction(u,phi,bc,costOldPhi);
                        costOldPhi = costPhi;
                        obj.printCost('iterPhi',iterPhi,costPhi,ePhi);
                        iterPhi = iterPhi + 1;

                        costFun(1,end+1) = costPhi;
                        costFun(2,end) = 1;
                        iter = iter+1;
                        obj.monitor.update(iter,{[],[],[],[],[costFun(1,end)],[]});
                    end
                    if iterPhi > iterPhiMax
                        iterPhiMax = iterPhi;
                    end
                    

                    [eStag, costStag] = obj.computeErrorCostFunction(u,phi,bc,costOldStag);
                    costOldStag = costStag;
                    obj.printCost('iterStag',iterStag,costStag,eStag);
                    iterStag = iterStag + 1;

                    costFun(1,end+1) = costStag;
                    costFun(2,end) = 2;
                    iter = iter+1;
                    obj.monitor.update(iter,{[],[],[],[],[costFun(1,end)],[]});
                end
                uOld = u;
                phiOld = phi;


                %%% SAVE DATA + MONITORING %%%%%
                s.step = i;
                s.numIterU = iterUMax-1;
                s.numIterP = iterPhiMax-1;
                s.numIterStag = iterStag-1;
                s.u = u; s.phi = phi; s.F = F;
                s.bc = bc;
                s.damageField = phi;
                s.cost = costFun;
                output = obj.saveData(output,s);
                
                totE = obj.functional.computeCostFunction(u,phi,bc);
                totF = obj.computeTotalReaction(F);
                displ = obj.boundaryConditions.bcValues(i);
                obj.monitor.update(i,{[totF;displ],[max(phi.fValues);displ],[phi.fValues],...
                                    [iterStag-1],[],[totE]});

            end
        end

    end


    methods (Access = private)
        %% %%%%%%%%%%%%%%%%%%%%%% INITIALIZATION %%%%%%%%%%%%%%%%%%%%%%%% %%
        function init(obj,cParams)
            obj.mesh                = cParams.mesh;
            obj.initPhi             = cParams.initPhi;
            obj.boundaryConditions  = cParams.boundaryConditions;
            obj.functional          = cParams.functional;
            obj.shallPrint          = cParams.monitoring.print;
            obj.setMonitoring(cParams)
        end

        function setMonitoring(obj,cParams)
            s.shallDisplay = cParams.monitoring.set;
            s.type = cParams.monitoring.type;
            s.mesh = obj.mesh;
            obj.monitor = PhaseFieldMonitoring.initialize(s);
        end

        function u = updateInitialDisplacement(obj,bc,uOld)
            restrictedDofs = bc.dirichlet_dofs;
            if isempty(restrictedDofs)
                u = uOld;
            else
                dirich = bc.dirichletFun;
                uVec = reshape(uOld.fValues',[uOld.nDofs 1]);
                dirichVec = reshape(dirich.fValues',[dirich.nDofs 1]);

                uVec(restrictedDofs) = dirichVec(restrictedDofs);
                u = reshape(uVec,[flip(size(uOld.fValues))])';
            end
        end

        %% %%%%%%%%%%%%%%%%%%%%%% ELASTIC EQUATION %%%%%%%%%%%%%%%%%%%%%%%% %%

        function uOut = computeDisplacement(obj,LHSfull, RHSfull,uIn,bc)
            [LHS,RHS] = fullToReduced(obj,LHSfull,RHSfull,bc);
            if ~isempty(LHS)
                uInVec = reshape(uIn.fValues',[uIn.nDofs 1]);
                uOutVec = uInVec;

                uInFree = uInVec(bc.free_dofs);
                uOutFree = obj.updateWithNewton(LHS,RHS,uInFree);
                uOutVec(bc.free_dofs) = uOutFree;
                uOut = reshape(uOutVec,[flip(size(uIn.fValues))])';
            else
                uOut = uIn.fValues;
            end
        end

        function [LHS,RHS] = fullToReduced(obj,LHS,RHS,bc)
            free_dofs = bc.free_dofs;
            LHS = LHS(free_dofs, free_dofs);
            RHS = RHS(free_dofs);
        end

        function xNew = updateWithNewton(obj,LHS,RHS,x)
            deltaX = -LHS\RHS;
         %   deltaX = -obj.tau.*RHS;
            xNew = x + deltaX; 
        end

        function F = computeForceVector(obj,LHS,u)
            uVec = reshape(u.fValues',[u.nDofs 1]);
            F = LHS*uVec;
        end

        %% %%%%%%%%%%%%%%%% PHASE-FIELD EQUATION %%%%%%%%%%%%%%%%%%%%%%%% %%

        function xP = projectInLowerAndUpperBound(obj,x,xLB,xUB)
            xP = min(max(xLB, x),xUB);
        end

        function xNew = updatePhi(obj,LHS,RHS,x)
            deltaX = -LHS\RHS;
          % deltaX = -obj.tau.*RHS;
            xNew = x + deltaX; 
        end

        %% %%%%%%%%%%%%%%%%%% AUXILIARY METHODS %%%%%%%%%%%%%%% %%

        function totReact = computeTotalReaction(obj,F)
            UpSide  = max(obj.mesh.coord(:,2));
            isInUp = abs(obj.mesh.coord(:,2)-UpSide)< 1e-12;
            nodes = 1:obj.mesh.nnodes;
            totReact = sum(F(2*nodes(isInUp)-1));

            % DownSide  = min(obj.mesh.coord(:,2));
            % isInDown = abs(obj.mesh.coord(:,2)-DownSide)< 1e-12;
            % nodes = 1:obj.mesh.nnodes;
            % totReact = sum(F(2*nodes(isInDown)));
        end

        function [e, cost] = computeErrorCostFunction(obj,u,phi,bc,costOld)
            cost = obj.functional.computeCostFunction(u,phi,bc);
            e = abs(cost - costOld);
        end

        function printCost(obj,name,iter,cost,e)
            if obj.shallPrint == true
                X = sprintf('%s:%d / cost: %.8e  (diff:%.8e) \n',name,iter,cost,e);
                fprintf(X);
            end
        end

        %% %%%%%%%%%%%%%%%%%% SAVE + PLOTS %%%%%%%%%%%%%%% %%

        function data = saveData(obj,data,cParams)
            step = cParams.step;
            fExt = cParams.bc.pointloadFun;

            data.reaction(step) = obj.computeTotalReaction(cParams.F);
            data.displacement(step) = obj.boundaryConditions.bcValues(step);
            data.damage.maxValue(step) = max(cParams.phi.fValues);
            data.damage.field = cParams.phi;

            data.energy.extWork(step) = -obj.functional.computeExternalWork(cParams.u,fExt);
            data.energy.intE(step) = obj.functional.computeInternalEnergy(cParams.u,cParams.phi);
            data.energy.localDis(step) = obj.functional.computeDissipationEnergy(cParams.phi);
            data.energy.regDis(step) = obj.functional.computeRegularisationEnergy(cParams.phi);
            data.iter.u(step) = cParams.numIterU;
            data.iter.phi(step) = cParams.numIterP;
            data.iter.stag(step) = cParams.numIterStag;
            data.cost = cParams.cost;
        end

    end

end