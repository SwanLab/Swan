classdef PhaseFieldComputer < handle

    properties (Constant, Access = public)
        tolErrU = 1e-12;
        tolErrPhi = 1e-12;
        tolErrStag = 1e-12;
        tau = 1e2;
    end

    properties (Access = private)
        steps
        bcVal 

        mesh
        BC
        
        materialPhaseField
        dissipationInterpolation

        fem
        phi
        u
        forceVector
    end

    properties (Access = private)
        functional
        data

        constant
        l0

        forceMat
        reactionMat
        displacementMat
        damageMat
        damageTheory
        stressStrainMat
        energyMat
        iterMat
        costFun

        IntegralCost
        IntegralEnergy
        IntegralDissipation
    end

    methods (Access = public)

        function obj = PhaseFieldComputer(cParams)
            obj.init(cParams)
            obj.createFunctionals()
            
            u = LagrangianFunction.create(obj.mesh,2,'P1');
            phi = LagrangianFunction.create(obj.mesh,1,'P1');
            uOld = u;
            phiOld = phi;

            obj.costFun = null(2,1);

            for i = 1:obj.steps
                disp([newline '%%%%%%%%%% STEP ',num2str(i),' %%%%%%%%%%%'])
                obj.createBoundaryConditions(obj.bcVal(i));
                u.fValues = obj.updateInitialDisplacement(uOld);

                obj.costFun(1,end+1) = obj.computeCostFunction(u,phi);
                obj.costFun(2,end) = 2;

                eStag = 1; iterStag = 1; costOldStag = 0;
                iterUMax = 1; iterPhiMax = 1;
                while (eStag > obj.tolErrStag) && (iterStag < 100)

                    eU = 1; iterU = 1; costOldU = 0;
                    while (eU > obj.tolErrU) && (iterU < 100)
                        LHS = obj.computeElasticLHS(u,phi);
                        RHS = obj.computeElasticResidual(u,phi);
                        u.fValues = obj.computeDisplacement(LHS,RHS,u);
                        F = obj.computeForces(LHS,u);

                        [eU, costU] = obj.computeErrorCostFunction(u,phi,costOldU);
                        costOldU = costU;
                        disp(['iterU: ',num2str(iterU),' res: ',num2str(eU)])

                        iterU = iterU + 1;
                    end
                    if iterU > iterUMax
                        iterUMax = iterU;
                    end

                    ePhi = 1;  iterPhi = 1; costOldPhi = 0;
                    while (ePhi > obj.tolErrPhi) && (iterPhi < 100)
                        LHS = obj.computePhaseFieldLHS(u,phi);
                        RHS = obj.computePhaseFieldResidual(u,phi);
                        phiNew = obj.updateWithNewton(LHS, RHS, phi.fValues);
                        phi.fValues = obj.projectInLowerAndUpperBound(phiNew,phiOld.fValues,1);

                        [ePhi, costPhi] = obj.computeErrorCostFunction(u,phi,costOldPhi);
                        costOldPhi = costPhi;
                        disp(['iterPhi: ',num2str(iterPhi),' res: ',num2str(ePhi)])

                        iterPhi = iterPhi + 1;
                    end
                    if iterPhi > iterPhiMax
                        iterPhiMax = iterPhi;
                    end

                    [eStag, costStag] = obj.computeErrorCostFunction(u,phi,costOldStag);
                    costOldStag = costStag;
                    disp(['iterStag: ',num2str(iterStag),' res: ',num2str(eStag)])

                    iterStag = iterStag + 1;
                end
                uOld = u;
                phiOld = phi;

                s.step = i;
                s.numIterU = iterUMax-1;
                s.numIterP = iterPhiMax-1;
                s.numIterStag = iterStag-1;
                s.u = u; s.phi = phi; s.F = F;
                obj.saveData(s);

                figure(400)
                hold on
                phi.plot;
                title(['Damage at step ',num2str(i)])
                colorbar
                clim([0 1])

                %obj.printPlots(phi,i);
                %exportgraphics(gcf,'Lshape.gif','Append',true);
            end
            obj.printPlots(phi,i);
        end

    end


    methods (Access = private)
        %% %%%%%%%%%%%%%%%%%%%%%% INITIALIZATION %%%%%%%%%%%%%%%%%%%%%%%% %%
        function init(obj,cParams)
            obj.bcVal = cParams.bcVal;
            obj.steps = length(cParams.bcVal);

            obj.mesh                     = cParams.mesh;
            obj.materialPhaseField       = cParams.materialPhaseField;
            obj.dissipationInterpolation = cParams.dissipationPhaseField;
            obj.l0 = cParams.l0;
            obj.constant = cParams.constant;

            obj.data.mat.reaction = zeros(1,obj.steps);
            obj.data.mat.displacementTip = zeros(1,obj.steps);
            obj.data.mat.maxDamage = zeros(1,obj.steps);
            obj.data.mat.stress = zeros(1,obj.steps);
            obj.data.mat.strain = zeros(1,obj.steps);
            obj.data.mat.energies = zeros(4,obj.steps);
            obj.data.mat.iterations = zeros(3,obj.steps);
        end

        function createBoundaryConditions(obj,prescribedVal)
            s.mesh = obj.mesh;
            bc = BoundaryContionsForPhaseFieldCreator(s);
            bcSet = bc.create(prescribedVal);
            obj.BC = bcSet;
        end

        function createFunctionals(obj)
            s.mesh = obj.mesh;
            s.materialPhaseField = obj.materialPhaseField;
            s.dissipationInterpolation = obj.dissipationInterpolation;
            s.constant = obj.constant;
            s.l0 = obj.l0;

            obj.functional.energy         = ShFunc_InternalEnergy(s);
            obj.functional.localDamage    = ShFunc_LocalDamage(s);
            obj.functional.nonLocalDamage = ShFunc_NonLocalDamage(s);
            obj.functional.extWork        = ShFunc_ExternalWork(s);
        end

        function u = updateInitialDisplacement(obj,u)
            restrictedDofs = obj.BC.dirichlet_dofs;
            dirich = obj.BC.dirichletFun;
            uVec = reshape(u.fValues',[u.nDofs 1]);
            dirichVec = reshape(dirich.fValues',[dirich.nDofs 1]);

            uVec(restrictedDofs) = dirichVec(restrictedDofs);
            u = reshape(uVec,[flip(size(u.fValues))])';
        end

        %% %%%%%%%%%%%%%%%%%%%%%% ELASTIC EQUATION %%%%%%%%%%%%%%%%%%%%%%%% %%

        function RHS = computeElasticResidual(obj,u,phi)
            fExt = obj.BC.pointloadFun;
            [Fint,~] = obj.functional.energy.computeGradient(u,phi,2);
            Fext = obj.functional.extWork.computeGradient(u,fExt,1);
            RHS = Fint + Fext;
        end

        function LHS = computeElasticLHS(obj,u,phi)
            [LHS,~] = obj.functional.energy.computeHessian(u,phi,2);
        end

        function uOut = computeDisplacement(obj,LHSfull, RHSfull,uIn)
            [LHS,RHS] = fullToReduced(obj,LHSfull,RHSfull);
            if ~isempty(LHS)
                uInVec = reshape(uIn.fValues',[uIn.nDofs 1]);
                uOutVec = uInVec;

                uInFree = uInVec(obj.BC.free_dofs);
                uOutFree = obj.updateWithNewton(LHS,RHS,uInFree);
                uOutVec(obj.BC.free_dofs) = uOutFree;
                uOut = reshape(uOutVec,[flip(size(uIn.fValues))])';
            else
                uOut = uIn.fValues;
            end
        end

        function [LHS,RHS] = fullToReduced(obj,LHS,RHS)
            free_dofs = obj.BC.free_dofs;
            LHS = LHS(free_dofs, free_dofs);
            RHS = RHS(free_dofs);
        end

        function F = computeForces(obj,LHS,u)
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
        end

        function totVal = computeTotalDissipationLocal(obj,phi)
            totVal = obj.functional.localDamage.computeFunction(phi,2);
        end

        function totVal = computeTotalRegularizationTerm(obj,phi)
            totVal = obj.functional.nonLocalDamage.computeFunction(phi,2);
        end

        function totVal = computeTotalExternalWork(obj,u)
            fExt = obj.BC.pointloadFun;
            totVal = obj.functional.extWork.computeFunction(u,fExt,2);
        end

        %% %%%%%%%%%%%%%%%%%% AUXILIARY METHODS %%%%%%%%%%%%%%% %%

        function xNew = updateWithNewton(obj,LHS,RHS,x)
            deltaX = -LHS\RHS;
            xNew = x + deltaX; 
        end

        function totReact = computeTotalReaction(obj,F)
            % UpSide  = max(obj.mesh.coord(:,2));
            % isInUp = abs(obj.mesh.coord(:,2)-UpSide)< 1e-12;
            % nodes = 1:obj.mesh.nnodes;
            % totReact = -sum(F(2*nodes(isInUp)));

            DownSide  = min(obj.mesh.coord(:,2));
            isInDown = abs(obj.mesh.coord(:,2)-DownSide)< 1e-12;
            nodes = 1:obj.mesh.nnodes;
            totReact = sum(F(2*nodes(isInDown)));
        end

        function [e, cost] = computeErrorCostFunction(obj,u,phi,costOld)
            cost = computeCostFunction(obj,u,phi);
            e = abs(cost - costOld);
        end

        function Etot = computeCostFunction(obj,u,phi)
            Wext = obj.computeTotalExternalWork(u);
            Eint = obj.computeTotalInternalEnergy(u,phi);
            Edis = obj.computeTotalDissipationLocal(phi);
            Ereg = obj.computeTotalRegularizationTerm(phi);
            Etot = Eint+Edis+Ereg-Wext;
        end

        function computeEnergyIntegrals(obj,step)
            fMat = [0 obj.reactionMat];
            dMat = [0 obj.displacementMat];

            obj.IntegralCost(step) = trapz(dMat(1:step+1),fMat(1:step+1));
            obj.IntegralEnergy(step) = 0.5*(dMat(step+1)*fMat(step+1));
            obj.IntegralDissipation = obj.IntegralCost - obj.IntegralEnergy;
        end

        %% %%%%%%%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%%% %%

        function saveData(obj,cParams)
            step = cParams.step;

            %obj.forceMat(step) = obj.fem.stressFun.fValues(2,1,1);  %% ONLY ONE ELEMENT
            obj.reactionMat(step) = -obj.computeTotalReaction(cParams.F);
            
            obj.displacementMat(step) = obj.bcVal(step);
            obj.damageMat(step) = max(cParams.phi.fValues);
         
            % quad = Quadrature.set(obj.mesh.type);
            % quad.computeQuadrature('QUADRATIC');
            % obj.materialPhaseField.computeMatIso(quad);
            % C = obj.materialPhaseField.material.C(2,2,1,5);
            % Gc = obj.materialPhaseField.Gc;
            % e = obj.displacementMat(step);
            % obj.damageTheory(step) = (C*e^2)/((Gc/obj.l0)+(C*e^2));

            obj.energyMat(1,step) = -obj.computeTotalExternalWork(cParams.u);
            obj.energyMat(2,step) = obj.computeTotalInternalEnergy(cParams.u,cParams.phi);
            obj.energyMat(3,step) = obj.computeTotalDissipationLocal(cParams.phi);
            obj.energyMat(4,step) = obj.computeTotalRegularizationTerm(cParams.phi);
            obj.energyMat(5,step) = sum(obj.energyMat(2,step) + obj.energyMat(1,step));

            obj.iterMat(1,step) = cParams.numIterU;              
            obj.iterMat(2,step) = cParams.numIterP;
            obj.iterMat(3,step) = cParams.numIterStag;

            obj.computeEnergyIntegrals(step)
        end

        function printPlots(obj,phi,step)
            figure(400)
            hold on
            phi.plot;
            title(['Damage at step ',num2str(step)])
            colorbar
            clim([0 1])

            % figure(100)
            % plot(obj.displacementMat(1:step),obj.forceMat(1:step))
            % title('Force-displacement diagram (Internal force)')
            % grid on
            % xlabel('Displacement [mm]')
            % ylabel('Force [kN]')

            figure(101)
            plot(obj.displacementMat(1:step),obj.reactionMat(1:step))
            title('Force-displacement diagram (Reaction force)')
            grid on
            xlabel('Displacement [mm]')
            ylabel('Force [kN]')

            figure(102)
            hold on
            plot(obj.displacementMat(1:step),obj.damageMat(1:step),'Color',"#0072BD")
            %plot(obj.displacementMat(1:step),obj.damageTheory(1:step),'--','Color',"#D95319")
            title('Damage-displacement diagram')
            grid on
            % legend('Algorithm result', ...
            %        'Theoretical result')
            xlabel('Displacement [mm]')
            ylabel('Damage [-]')
            ylim([0 1]);

            figure(200)
            hold on
            plot(obj.energyMat(5,1:step),'Color',"#0072BD")
            plot(obj.energyMat(2,1:step),'Color',"#D95319")
            plot(obj.energyMat(3,1:step),'Color',"#EDB120")
            plot(obj.energyMat(4,1:step),'Color',"#7E2F8E")
            plot(obj.energyMat(1,1:step))
            title('Energy values at each step (Equation terms)')
            legend('Total energy', ...
                   'Internal Energy', ...
                   'Local surface energy', ...
                   'Non-local surface energy', ...
                   'External Work')
            xlabel('Step [-]')
            ylabel('Energy [J]')

            figure(300)
            hold on
            plot(obj.iterMat(1,1:step),'xb')
            plot(obj.iterMat(2,1:step),'r')
            plot(obj.iterMat(3,1:step),'k')
            title('Iterations needed')
            legend('U','phi')
            xlabel('Step [-]')
            ylabel('Iterations')


            % figure(500)
            % hold on
            % for n = 2:size(obj.costFun,2)
            %     if obj.costFun(2,n) == 0
            %         plot([n-1, n],[obj.costFun(1,n-1), obj.costFun(1,n)],'b')
            %     elseif obj.costFun(2,n) == 1
            %         plot([n-1, n],[obj.costFun(1,n-1), obj.costFun(1,n)],'r')
            %     elseif obj.costFun(2,n) == 2
            %         plot([n-1, n],[obj.costFun(1,n-1), obj.costFun(1,n)],'k')
            %     end
            % end
            % title('Cost Function')
            % xlabel('Iteration [-]')
            % ylabel('Energy [J]')

            figure(600)
            hold on
            plot(obj.IntegralCost(1:step))
            plot(obj.IntegralEnergy(1:step))
            plot(obj.IntegralDissipation(1:step))
            title('Energy values at each step (Force-displacement)')
            legend('Total energy', ...
                   'Internal energy', ...
                   'Dissipated energy')
            xlabel('Step [-]')
            ylabel('Energy [J]')

        end

        function printResults(obj,step)
            obj.fem.print(['PhaseFieldFEM_Step',step],'GiD');
            obj.phi.print(['PhaseFieldDamage_Step',step],'GiD');
        end

    end

end
