classdef PhaseFieldComputer < handle

    properties (Access = private)
        mesh
        boundaryConditions
        initialGuess
        functional
        tol
        solverType
    end

    properties (Access = private)
        shallPrint
        monitor
        stop
    end

    methods (Access = public)

        function obj = PhaseFieldComputer(cParams)
            obj.init(cParams)
        end

        function outputData = compute(obj)
            u   = obj.initialGuess.u;
            phi = obj.initialGuess.phi;
            uOld = u;
            phiOld = phi;

            iter = 0;
            phiIter = 0;
            costFun = 0;
            tau = 150;
            tauArray = [];

            maxSteps = length(obj.boundaryConditions.bcValues);
            obj.stop.noFullyBroken = true; obj.stop.triggered = false;
            obj.stop.maxF = 0;             obj.stop.lastStep = maxSteps;
            i = 1;
            while(i<=maxSteps) && (obj.stop.noFullyBroken)
                obj.monitor.printStep(i,maxSteps)
                bc = obj.boundaryConditions.nextStep();
                u.setFValues(obj.updateInitialDisplacement(bc,uOld));
                % updateLowerBound()
                eStag = 1; iterStag = 1; costOldStag = 0;
                iterUMax = 1; iterPhiMax = 1;
                while (abs(eStag) > obj.tol.stag) && (iterStag < 300)

                    cParams.functional = obj.functional;
                    cParams.tolerance = obj.tol.u;
                    cParams.maxIter  = iterUMax;
                    cParams.monitor  = obj.monitor;
                    A = PhaseFieldDisplacementUpdater(cParams);
                    [u,F,costFun,iterU] = A.update(u,phi,bc,costFun);

                    if iterU > iterUMax
                        iterUMax = iterU;
                    end

                    ePhi = -1;  iterPhi = 1; costOldPhi = costFun(end);
                    while (abs(ePhi) > obj.tol.phi) && (iterPhi < 300)
                        LHS = obj.functional.computePhaseFieldLHS(u,phi);
                        RHS = obj.functional.computePhaseFieldRHS(u,phi);

                        if obj.solverType == "Gradient"
                            % ADAPTIVE LINE-SEARCH GRADIENT
                            phiNew = obj.updateWithGradient(RHS, phi.fValues,tau);
                            phiProposed = phi.copy();
                            phiProposed.setFValues(obj.projectInLowerAndUpperBound(phiNew,phiOld.fValues,1));
                            [ePhi, costPhi] = obj.computeErrorCostFunctional(u,phiProposed,bc,costOldPhi);

                            obj.monitor.printCost('iterPhi',iterPhi,costPhi,ePhi);
                            iterPhi = iterPhi + 1;

                            if ePhi > 0
                                tau = tau/2;
                            else
                                phiIter = phiIter + 1;
                                obj.monitor.update(phiIter,{[],[],[],[],[],[],[tau]})
                                tauArray(end+1) = tau;
                                if tau<=1e10
                                    tau = 10*tau;
                                end
                                phi = phiProposed;
                                costOldPhi = costPhi;
                                costFun(end+1) = costPhi;
                                iter = iter+1;
                                obj.monitor.update(iter,{[],[],[],[],[costFun(end)],[],[]});
                            end



                        elseif obj.solverType == "Newton"
                            %NEWTON METHOD
                            phiNew = obj.updateWithNewton(LHS,RHS,phi.fValues);
                            phi.setFValues (obj.projectInLowerAndUpperBound(phiNew,phiOld.fValues,1));
                            [ePhi, costPhi] = obj.computeErrorCostFunctional(u,phi,bc,costOldPhi);

                            costOldPhi = costPhi;
                            obj.monitor.printCost('iterPhi',iterPhi,costPhi,ePhi);
                            iterPhi = iterPhi + 1;

                            costFun(end+1) = costPhi;
                            iter = iter+1;
                            obj.monitor.update(iter,{[],[],[],[],[costFun(end)],[],[]});
                            obj.monitor.refresh();
                        end
                    end

                    if iterPhi > iterPhiMax
                        iterPhiMax = iterPhi;
                    end


                    [eStag, costStag] = obj.computeErrorCostFunctional(u,phi,bc,costOldStag);
                    costOldStag = costStag;
                    obj.monitor.printCost('iterStag',iterStag,costStag,eStag);
                    iterStag = iterStag + 1;

                    costFun(end+1) = costStag;
                    iter = iter+1;
                    obj.monitor.update(iter,{[],[],[],[],[costFun(end)],[],[]});
                    obj.monitor.refresh();
                end
                uOld = u;
                phiOld = phi;

                %%% SAVE DATA + MONITORING %%%%% KEEP CLEANING
                totE = obj.functional.computeCostFunctional(u,phi,bc);
                totF = obj.computeTotalReaction(F);

                s.force = totF;
                s.bcVal = obj.boundaryConditions.bcValues(i);
                fExt = bc.pointloadFun;
                s.energy = obj.functional.computeEnergiesFunctional(u,phi,fExt);
                s.numIterU = iterUMax-1;
                s.numIterP = iterPhiMax-1;
                s.numIterStag = iterStag-1;
                s.u = u; s.phi = phi; s.F = F;
                s.bc = bc;
                s.damageField = phi;
                s.cost = costFun;
                s.tauArray = tauArray;
                obj.monitor.saveData(i,s);


                displ = obj.boundaryConditions.bcValues(i);
                obj.monitor.update(i,{[totF;displ],[max(phi.fValues);displ],[phi.fValues],...
                    [iterStag-1],[],[totE;displ],[]});
                obj.monitor.refresh();



                obj.checkStopCode(i,totF);
                i = i + 1;



            end
            outputData = obj.monitor.data;
        end

    end

    methods (Access = private)
        %% %%%%%%%%%%%%%%%%%%%%%% INITIALIZATION %%%%%%%%%%%%%%%%%%%%%%%% %%
        function init(obj,cParams)
            obj.mesh               = cParams.mesh;
            obj.initialGuess       = cParams.initialGuess;
            obj.boundaryConditions = cParams.boundaryConditions;
            obj.functional         = cParams.functional;
            obj.shallPrint         = cParams.monitoring.print;
            obj.tol                = cParams.tolerance;
            obj.solverType         = cParams.solverType;
            obj.setMonitoring(cParams)
        end

        function setMonitoring(obj,cParams)
            s.shallDisplay = cParams.monitoring.set;
            s.shallPrint   = cParams.monitoring.print;
            s.type = cParams.monitoring.type;
            s.fun = obj.initialGuess.phi;
            obj.monitor = PhaseFieldMonitoring(s);
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

        function checkStopCode(obj,i,totF)
            if totF > obj.stop.maxF
                obj.stop.maxF = totF;
            elseif i>5 && totF<0.01*obj.stop.maxF && ~obj.stop.triggered
                obj.stop.lastStep = i;
                obj.stop.triggered = true;
            end

            if i==obj.stop.lastStep+10
                obj.stop.noFullyBroken = false;
            end
        end


        %% %%%%%%%%%%%%%%%% PHASE-FIELD EQUATION %%%%%%%%%%%%%%%%%%%%%%%% %%

        function xP = projectInLowerAndUpperBound(obj,x,xLB,xUB)
            xP = min(max(xLB, x),xUB);
        end

        function xNew = updateWithNewton(obj,LHS,RHS,x)
            deltaX = -LHS\RHS;
            xNew = x + deltaX;
        end

        function xNew = updateWithGradient(obj,RHS,x,tau)
            deltaX = -tau.*RHS;
            xNew = x + deltaX;
        end

        %% %%%%%%%%%%%%%%%%%% AUXILIARY METHODS %%%%%%%%%%%%%%% %%

        function totReact = computeTotalReaction(obj,F)
            UpSide  = max(obj.mesh.coord(:,2));
            isInUp = abs(obj.mesh.coord(:,2)-UpSide)< 1e-12;
            nodes = 1:obj.mesh.nnodes;
            ReactX = sum(F(2*nodes(isInUp)-1));
            ReactY = sum(F(2*nodes(isInUp)));
            totReact = sqrt(ReactX^2+ReactY^2);

            % DownSide  = min(obj.mesh.coord(:,2));
            % isInDown = abs(obj.mesh.coord(:,2)-DownSide)< 1e-12;
            % nodes = 1:obj.mesh.nnodes;
            % totReact = -sum(F(2*nodes(isInDown)));

            % isInTip = (abs(obj.mesh.coord(:,2)-(max(obj.mesh.coord(:,2))+min(obj.mesh.coord(:,2)))/2) < 1e-12) & (abs(obj.mesh.coord(:,1)-max(obj.mesh.coord(:,1))) < 30);
            % nodes = 1:obj.mesh.nnodes;
            % totReact = sum(F(2*nodes(isInTip)));
        end

        function [e, cost] = computeErrorCostFunctional(obj,u,phi,bc,costOld)
            cost = obj.functional.computeCostFunctional(u,phi,bc);
            e = cost - costOld;
        end





    end
end
