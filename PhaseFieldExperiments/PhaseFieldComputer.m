classdef PhaseFieldComputer < handle

    properties (Constant, Access = public)
        tolErrU = 1e-6;
        tolErrPhi = 1e-8;
        tolErrStag = 1e-8;
        tau = 100*1e2;
    end

    properties (Access = private)
        mesh
        boundaryConditions
        initPhi
        material
        dissipation

    end

    properties (Access = private)
        monitor
        functional
        l0

        costFun
        data
    end

    methods (Access = public)

        function obj = PhaseFieldComputer(cParams)
            obj.init(cParams)
            obj.createFunctionals()
            obj.setMonitoring(cParams)
        end

        function output = compute(obj)
            u = LagrangianFunction.create(obj.mesh,2,'P1');
            phi = obj.initPhi;
            uOld = u;
            phiOld = phi;

            obj.costFun = null(2,1);
            output = [];
            %costStag = 0;
            maxSteps = length(obj.boundaryConditions.bcValues);
            for i = 1:maxSteps
                disp([newline '%%%%%%%%%% STEP ',num2str(i),' %%%%%%%%%%%'])
                bc = obj.boundaryConditions.nextStep();
                u.fValues = obj.updateInitialDisplacement(bc,uOld);

                % obj.costFun(1,end+1) = costStag;
                % obj.costFun(2,end) = 2;

                eStag = 1; iterStag = 1; costOldStag = 0;
                iterUMax = 1; iterPhiMax = 1;
                while (eStag > obj.tolErrStag) && (iterStag < 300)

                    eU = 1; iterU = 1; costOldU = 0;
                    while (eU > obj.tolErrU) && (iterU < 100)
                        LHS = obj.computeElasticLHS(u,phi);
                        RHS = obj.computeElasticResidual(u,phi,bc);
                        u.fValues = obj.computeDisplacement(LHS,RHS,u,bc);
                        F = obj.computeForceVector(LHS,u);

                        [eU, costU] = obj.computeErrorCostFunction(u,phi,bc,costOldU);
                        costOldU = costU;
                        disp(['iterU: ',num2str(iterU),' res: ',num2str(eU)])

                        iterU = iterU + 1;
                    end
                    if iterU > iterUMax
                        iterUMax = iterU;
                    end
                   
                  % obj.costFun(1,end+1) = costU;
                  % obj.costFun(2,end) = 0;  
                  % obj.printPlots()


                    ePhi = 1;  iterPhi = 1; costOldPhi = 0;
                    while (ePhi > obj.tolErrPhi) && (iterPhi < 100)
                        LHS = obj.computePhaseFieldLHS(u,phi);
                        RHS = obj.computePhaseFieldResidual(u,phi);
                        phiNew = obj.updatePhi(LHS, RHS, phi.fValues);
                        phi.fValues = obj.projectInLowerAndUpperBound(phiNew,phiOld.fValues,1);

                        [ePhi, costPhi] = obj.computeErrorCostFunction(u,phi,bc,costOldPhi);
                        costOldPhi = costPhi;
                        disp(['iterPhi: ',num2str(iterPhi),' res: ',num2str(ePhi)])

                        iterPhi = iterPhi + 1;
                    end
                    if iterPhi > iterPhiMax
                        iterPhiMax = iterPhi;
                    end
                    
                    % obj.costFun(1,end+1) = costPhi;
                    % obj.costFun(2,end) = 1; 
                    % obj.printPlots()


                    [eStag, costStag] = obj.computeErrorCostFunction(u,phi,bc,costOldStag);
                    costOldStag = costStag;
                    disp(['iterStag: ',num2str(iterStag),' res: ',num2str(eStag)])

                    iterStag = iterStag + 1;
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
                if i==maxSteps
                    s.damageField = phi;
                end
                output = obj.saveData(output,s);
                
                
                TotalForce = obj.computeTotalReaction(F);
                Displacement = obj.boundaryConditions.bcValues(i);
                obj.monitor.update(i,{[TotalForce;Displacement],[max(phi.fValues);Displacement],iterStag-1,phi.fValues});

            end
        end

    end


    methods (Access = private)
        %% %%%%%%%%%%%%%%%%%%%%%% INITIALIZATION %%%%%%%%%%%%%%%%%%%%%%%% %%
        function init(obj,cParams)
            obj.mesh                = cParams.mesh;
            obj.initPhi             = cParams.initPhi;
            obj.material            = cParams.material;
            obj.dissipation         = cParams.dissipation;
            obj.boundaryConditions  = cParams.boundaryConditions;
            obj.l0                  = cParams.l0;
        end

        function createFunctionals(obj)
            s.mesh = obj.mesh;
            s.material = obj.material;
            s.dissipation = obj.dissipation;
            s.l0 = obj.l0;

            obj.functional.energy         = ShFunc_InternalEnergySplit(s);
            obj.functional.energy2        = ShFunc_InternalEnergy(s);            
            obj.functional.localDamage    = ShFunc_LocalDamage(s);
            obj.functional.nonLocalDamage = ShFunc_NonLocalDamage(s);
            obj.functional.extWork        = ShFunc_ExternalWork(s);
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

        function RHS = computeElasticResidual(obj,u,phi,bc)
            fExt     = bc.pointloadFun;
            [Fint,~] = obj.functional.energy.computeGradient(u,phi,2);
            [Fint2,~] = obj.functional.energy2.computeGradient(u,phi,2);
            nFu = norm(Fint(:)-Fint2(:))/norm(Fint2(:))
            Fext     = obj.functional.extWork.computeGradient(u,fExt,1);
            RHS = Fint + Fext;
        end

        function LHS = computeElasticLHS(obj,u,phi)
            [LHS,~] = obj.functional.energy.computeHessian(u,phi,2);
            [LHS2,~] = obj.functional.energy2.computeHessian(u,phi,2);
            nKu = norm(LHS(:)-LHS2(:))/norm(LHS2(:))           
        end

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

        function F = computeForceVector(obj,LHS,u)
            uVec = reshape(u.fValues',[u.nDofs 1]);
            F = LHS*uVec;
        end

        %% %%%%%%%%%%%%%%%% PHASE-FIELD EQUATION %%%%%%%%%%%%%%%%%%%%%%%% %%

        function xP = projectInLowerAndUpperBound(obj,x,xLB,xUB)
            xP = min(max(xLB, x),xUB); 
        end

        function LHS = computePhaseFieldLHS(obj,u,phi)
            Mi = obj.createInternalEnergyMassMatrix(u,phi);
            Md = obj.createDissipationMassMatrix(phi);
            K  = obj.createStiffnessMatrix(phi);
            LHS = Mi + Md + K;
        end

        function res = computePhaseFieldResidual(obj,u,phi)
            Fi = obj.createInternalEnergyForceVector(u,phi);
            Fd = obj.createDissipationForceVector(phi);
            DF = obj.createForceDerivativeVector(phi);            
            res = Fi + Fd + DF;
        end
        
        
        %%% LHS %%%
        % Internal energy mass matrix
        function Mi = createInternalEnergyMassMatrix(obj,u,phi)
            [~,Hphiphi] = obj.functional.energy.computeHessian(u,phi,2);
            [~,Hphiphi2] = obj.functional.energy2.computeHessian(u,phi,2);            
            nKphi = norm(Hphiphi2(:) - Hphiphi(:))/norm(Hphiphi2(:))
            Mi = Hphiphi;
        end

        % Dissipation mass matrix
        function Md = createDissipationMassMatrix(obj,phi)
            Md = obj.functional.localDamage.computeHessian(phi,2);
        end

        % Stiffness matrix
        function K = createStiffnessMatrix(obj,phi)
            K = obj.functional.nonLocalDamage.computeHessian(phi,2);
        end

        %%% RHS %%%
        % Internal energy force vector
        function Fi = createInternalEnergyForceVector(obj,u,phi)
            [~,Jphi] = obj.functional.energy.computeGradient(u,phi,2);
            [~,Jphi2] = obj.functional.energy2.computeGradient(u,phi,2);            
            nFphi = norm(Jphi(:)-Jphi2(:))/norm(Jphi(:))
            Fi = Jphi;
        end
        
        % Dissipation force vector
        function Fd = createDissipationForceVector(obj,phi)
            Fd = obj.functional.localDamage.computeGradient(phi,2);     
        end
       
        % Force derivative vector
        function DF = createForceDerivativeVector(obj,phi)
            DF = obj.functional.nonLocalDamage.computeGradient(phi,2);
        end

        %% %%%%%%%%%%%%%%%% COMPUTE TOTAL ENERGIES %%%%%%%%%%%%%%%%%%%%%%%% %%
        function totVal = computeTotalInternalEnergy(obj,u,phi)
            totVal = obj.functional.energy.computeFunction(u,phi,2);
            totVal2 = obj.functional.energy2.computeFunction(u,phi,2);            
            nE = norm(totVal2-totVal)/norm(totVal2)
        end

        function totVal = computeTotalDissipationLocal(obj,phi)
            totVal = obj.functional.localDamage.computeFunction(phi,2);
        end

        function totVal = computeTotalRegularizationTerm(obj,phi)
            totVal = obj.functional.nonLocalDamage.computeFunction(phi,2);
        end

        function totVal = computeTotalExternalWork(obj,u,bc)
            fExt = bc.pointloadFun;
            totVal = obj.functional.extWork.computeFunction(u,fExt,2);
        end

        %% %%%%%%%%%%%%%%%%%% AUXILIARY METHODS %%%%%%%%%%%%%%% %%

        function xNew = updateWithNewton(obj,LHS,RHS,x)
            deltaX = -LHS\RHS;
         %   deltaX = -obj.tau.*RHS;
            xNew = x + deltaX; 
        end

        function xNew = updatePhi(obj,LHS,RHS,x)
            deltaX = -LHS\RHS;
          % deltaX = -obj.tau.*RHS;
            xNew = x + deltaX; 
        end
        

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
            cost = computeCostFunction(obj,u,phi,bc);
            e = abs(cost - costOld);
        end

        function Etot = computeCostFunction(obj,u,phi,bc)
            Wext = obj.computeTotalExternalWork(u,bc);
            Eint = obj.computeTotalInternalEnergy(u,phi);
            Edis = obj.computeTotalDissipationLocal(phi);
            Ereg = obj.computeTotalRegularizationTerm(phi);
            Etot = Eint+Edis+Ereg-Wext;
        end

        %% %%%%%%%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%%% %%

        function data = saveData(obj,data,cParams)
            step = cParams.step;
            
            data.reaction(step) = obj.computeTotalReaction(cParams.F);
            data.displacement(step) = obj.boundaryConditions.bcValues(step);
            data.damage.maxValue(step) = max(cParams.phi.fValues);
            data.damage.field = cParams.phi;
            % quad = Quadrature.set(obj.mesh.type);
            % quad.computeQuadrature('QUADRATIC');
            % obj.materialPhaseField.computeMatIso(quad);
            % C = obj.materialPhaseField.material.C(2,2,1,5);
            % Gc = obj.materialPhaseField.Gc;
            % e = obj.displacementMat(step);
            % obj.damageTheory(step) = (C*e^2)/((Gc/obj.l0)+(C*e^2));

            data.energy.extWork(step) = -obj.computeTotalExternalWork(cParams.u,cParams.bc);
            data.energy.intE(step) = obj.computeTotalInternalEnergy(cParams.u,cParams.phi);
            data.energy.localDis(step) = obj.computeTotalDissipationLocal(cParams.phi);
            data.energy.regDis(step) = obj.computeTotalRegularizationTerm(cParams.phi);
            data.iter.u(step) = cParams.numIterU;              
            data.iter.phi(step) = cParams.numIterP;
            data.iter.stag(step) = cParams.numIterStag;
            
        end

        function printPlots(obj,phi,step)


            figure(300)
            hold on
            for n = 2:size(obj.costFun,2)
                if obj.costFun(2,n) == 0
                    plot([n-1, n],[obj.costFun(1,n-1), obj.costFun(1,n)],'b')
                elseif obj.costFun(2,n) == 1
                    plot([n-1, n],[obj.costFun(1,n-1), obj.costFun(1,n)],'r')
                elseif obj.costFun(2,n) == 2
                    plot([n-1, n],[obj.costFun(1,n-1), obj.costFun(1,n)],'k')
                end
            end
            title('Cost Function')
            xlabel('Iteration [-]')
            ylabel('Energy [J]')
            hold off
        end

    end

end
