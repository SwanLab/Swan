classdef PhaseFieldComputer < handle

    properties (Constant, Access = public)
        tolErrU = 1e-12;
        tolErrPhi = 1e-12;
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
            uOld = u.fValues;
            phiOld = phi.fValues;

            obj.costFun = null(2,1);

            for i = 1:obj.steps
                disp([newline '%%%%%%%%%% STEP ',num2str(i),' %%%%%%%%%%%'])
                obj.createBoundaryConditions(obj.bcVal(i));
                u.fValues = obj.BC.dirichletFun.fValues;

                obj.costFun(1,end+1) = obj.computeCostFunction(u,phi);
                obj.costFun(2,end) = 2;

                eStag = 1; iterStag = 1;
                while (eStag > obj.tolErrU) && (iterStag < 100)

                    eU = 1;iterU = 1;
                    while (eU > obj.tolErrU) && (iterU < 100)
                        LHSfull = obj.computeElasticLHS(u,phi);
                        RHSfull = obj.computeElasticResidual(u,phi);
                        [LHS,RHS] = fullToReduced(obj,LHSfull,RHSfull);
                        U = reshape(u.fValues',[u.nDofs 1]);
                        uNew = obj.updateWithNewton(LHS,RHS,U(freeDofs));
                        F = LHSfull*U;
                        u = reshape(uNew,[flip(size(u.fValues))])';

                        obj.costFun(1,end+1) = obj.computeCostFunction(u,phi);
                        obj.costFun(2,end) = 0;
                        eU = abs(obj.costFun(1,end)-obj.costFun(1,end-1));
                        disp(['iterPhi: ',num2str(numIterP),' res: ',num2str(errorPhi)])
                    end

                    ePhi = 1;  iterPhi = 1;
                    while (ePhi > obj.tolErrPhi) && (iterPhi < 100)
                        LHS = obj.computePhaseFieldLHS(u,phi);
                        RHS = obj.computePhaseFieldResidual(u,phi);
                        phiNew = obj.updateWithNewton(LHS,RHS,phi.fValues);
                        phi = obj.projectInLowerAndUpperBound(phiNew,phiOld,1);

                        obj.costFun(1,end+1) = obj.computeCostFunction(u,phi);
                        obj.costFun(2,end) = 1;
                        ePhi = abs(obj.costFun(1,end)-obj.costFun(1,end-1));
                        disp(['iterPhi: ',num2str(numIterP),' res: ',num2str(errorPhi)])

                        numIterP = numIterP + 1;
                    end



                    eStag = obj.computeDisplacementError(uOld);
                    disp(['iterU: ',num2str(numIterU),' res: ',num2str(eStag(end))])

                    numIterU = numIterU + 1;
                    uOld.fValues = obj.u.fValues;
                end
                phiOld = obj.phi.fValues;

                s.step = i;
                s.numIterU = numIterU;
                s.numIterP = numIterP;
                obj.saveData(s);
            end
            obj.printPlots(i);
        end

    end


    methods (Access = private)
        %% %%%%%%%%%%%%%%%%%%%%%% INITIALIZATION %%%%%%%%%%%%%%%%%%%%%%%% %%
        function init(obj,cParams)
            obj.bcVal = cParams.bcVal;
            obj.steps = length(cParams.bcVal);

            obj.mesh                     = cParams.mesh;
            obj.phi               = cParams.initialPhaseField;
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
            obj.data.mat.iterations = zeros(2,obj.steps);
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

        
        function xNew = updateWithNewton(LHS,RHS,x)
            deltaX = -LHS\RHS;
            xNew = x + deltaX; 
        end
        %% %%%%%%%%%%%%%%%%%%%%%% ELASTIC EQUATION %%%%%%%%%%%%%%%%%%%%%%%% %%

        function RHS = computeElasticResidual(obj,u,phi)
            fExt = obj.BC.pointloadFun;
            [Fint,~] = obj.functional.energy.computeGradient(u,phi,'LINEAR');
            Fext = obj.functional.extWork.computeGradient(u,fExt,'LINEAR');
            RHS = Fint + Fext;
        end

        function LHS = computeElasticLHS(obj,u,phi)
            [LHS,~] = obj.functional.energy.computeHessian(u,phi,'LINEAR');
        end

        function UfV = computeDisplacement(obj,LHSfull, RHSfull)

        end

        function [LHS,RHS] = fullToReduced(obj,LHS,RHS)
            free_dofs = obj.BC.free_dofs;
            LHS = LHS(free_dofs, free_dofs);
            RHS = RHS(free_dofs);
        end

        function F = computeReactions(obj,LHS)

        end

        %% %%%%%%%%%%%%%%%% PHASE-FIELD EQUATION %%%%%%%%%%%%%%%%%%%%%%%% %%

        function xP = projectInLowerAndUpperBound(obj,x,xLB,xUB)
            xP = min(max(xLB, x),xUB); 
        end

        function LHS = computePhaseFieldLHS(obj)
            Mi = obj.createInternalEnergyMassMatrix();
            Md = obj.createDissipationMassMatrix();
            K  = obj.createStiffnessMatrix();
            LHS = Mi + Md + K;
        end

        function res = computePhaseFieldResidual(obj)
            Fi = obj.createInternalEnergyForceVector();
            Fd = obj.createDissipationForceVector();
            DF = obj.createForceDerivativeVector();            
            res = Fi + Fd + DF;
        end
        
        
        %%% LHS %%%
        % Internal energy mass matrix
        function Mi = createInternalEnergyMassMatrix(obj)
            [~,Hphiphi] = obj.functional.energy.computeHessian(obj.u,obj.phi,'QUADRATIC');
            Mi = Hphiphi;
        end

        % Dissipation mass matrix
        function Md = createDissipationMassMatrix(obj)
            Md = obj.functional.localDamage.computeHessian(obj.phi,'LINEAR');
        end

        % Stiffness matrix
        function K = createStiffnessMatrix(obj)
            K = obj.functional.nonLocalDamage.computeHessian(obj.phi,'LINEAR');
        end

        %%% RHS %%%
        % Internal energy force vector
        function Fi = createInternalEnergyForceVector(obj)
            [~,Jphi] = obj.functional.energy.computeGradient(obj.u,obj.phi,'QUADRATIC');
            Fi = Jphi;
        end
        
        % Dissipation force vector
        function Fd = createDissipationForceVector(obj)
            Fd = obj.functional.localDamage.computeGradient(obj.phi,'LINEAR');     
        end
       
        % Force derivative vector
        function DF = createForceDerivativeVector(obj)
            DF = obj.functional.nonLocalDamage.computeGradient(obj.phi,'LINEAR');
        end

        %% %%%%%%%%%%%%%%%% COMPUTE TOTAL ENERGIES %%%%%%%%%%%%%%%%%%%%%%%% %%
        function totVal = computeTotalInternalEnergy(obj)
            totVal = obj.functional.energy.computeFunction(obj.u,obj.phi,'QUADRATIC');
        end

        function totVal = computeTotalDissipationLocal(obj)
            totVal = obj.functional.localDamage.computeFunction(obj.phi,'QUADRATIC');
        end

        function totVal = computeTotalRegularizationTerm(obj)
            totVal = obj.functional.nonLocalDamage.computeFunction(obj.phi,'CONSTANT');
        end

        function totVal = computeTotalExternalWork(obj)
            fExt = obj.BC.pointloadFun;
            totVal = obj.functional.extWork.computeFunction(obj.u,fExt,'LINEAR');
        end

        %% %%%%%%%%%%%%%%%%%% AUXILIARY METHODS %%%%%%%%%%%%%%% %%

        function totReact = computeTotalReaction(obj,F)
            UpSide  = max(obj.mesh.coord(:,2));
            isInUp = abs(obj.mesh.coord(:,2)-UpSide)< 1e-12;
            nodes = 1:obj.mesh.nnodes;
            totReact = -sum(F(2*nodes(isInUp)));
        end

        function e = computeCostFunction(obj)
            eW = obj.computeTotalExternalWork();
            eT = obj.computeTotalInternalEnergy();
            eD = obj.computeTotalDissipationLocal();
            eR = obj.computeTotalRegularizationTerm();
            e = eW+eT+eD+eR;
        end

        function error = computeDisplacementError(obj,f)
            s.mesh = obj.mesh;
            s.fValues = obj.u.fValues - f.fValues;
            f = P1Function(s);
            int = Integrator.create('ScalarProduct',obj.mesh,'CONSTANT');
            error = sqrt(int.compute(f,f));
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
            obj.reactionMat(step) = -obj.computeTotalReaction();
            
            obj.displacementMat(step) = obj.bcVal(step);
            obj.damageMat(step) = max(obj.phi.fValues);
         
            % quad = Quadrature.set(obj.mesh.type);
            % quad.computeQuadrature('QUADRATIC');
            % obj.materialPhaseField.computeMatIso(quad);
            % C = obj.materialPhaseField.material.C(2,2,1,5);
            % Gc = obj.materialPhaseField.Gc;
            % e = obj.displacementMat(step);
            % obj.damageTheory(step) = (C*e^2)/((Gc/obj.l0)+(C*e^2));

            obj.energyMat(1,step) = -obj.computeTotalExternalWork();
            obj.energyMat(2,step) = obj.computeTotalInternalEnergy();
            obj.energyMat(3,step) = obj.computeTotalDissipationLocal();
            obj.energyMat(4,step) = obj.computeTotalRegularizationTerm();
            obj.energyMat(5,step) = sum(obj.energyMat(2,step) + obj.energyMat(1,step));

            obj.iterMat(1,step) = cParams.numIterU;              
            obj.iterMat(2,step) = cParams.numIterP;

            obj.computeEnergyIntegrals(step)
        end

        function printPlots(obj,step)
            figure(400)
            hold on
            obj.phi.plot;
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
            plot(obj.iterMat(1,1:step),'b')
            plot(obj.iterMat(2,1:step),'r')
            title('Iterations needed')
            legend('U','phi')
            xlabel('Step [-]')
            ylabel('Iterations')
           

            figure(500)
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
