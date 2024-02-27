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
        boundaryConditions
        
        materialPhaseField
        dissipationInterpolation

        fem
        phi
        deltaPhi
        u
        deltaU
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

            obj.u = LagrangianFunction.create(obj.mesh,2,'P1');
            uOld = LagrangianFunction.create(obj.mesh,2,'P1');
            phiOld = obj.phi.fValues;
            obj.costFun = null(2,1);

            for i = 1:obj.steps
                disp([newline '%%%%%%%%%% STEP ',num2str(i),' %%%%%%%%%%%'])
                obj.createBoundaryConditions(obj.bcVal(i));
                numIterU = 1;
                errorU = 1;
                while (errorU(end) > obj.tolErrU) && (numIterU < 100)
                    obj.computeFEM();
                    obj.costFun(1,end+1) = obj.computeCostFunction();
                    if numIterU == 1
                        obj.costFun(2,end) = 2;
                    else
                        obj.costFun(2,end) = 0;
                    end

                    numIterP = 1;
                    errorPhi = 1;
                    while (errorPhi > obj.tolErrPhi) && (numIterP < 100)
                        obj.computePhaseField();
                        obj.phi.fValues = obj.phi.fValues + obj.deltaPhi;
                        obj.phi.fValues = min(max(phiOld, obj.phi.fValues),1);

                        obj.costFun(1,end+1) = obj.computeCostFunction();
                        obj.costFun(2,end) = 1;

                        errorPhi = abs(obj.costFun(1,end)-obj.costFun(1,end-1));
                        disp(['iterPhi: ',num2str(numIterP),' res: ',num2str(errorPhi)])

                        numIterP = numIterP + 1;
                    end
                    errorU(end+1) = obj.computeDisplacementError(uOld);
                    disp(['iterU: ',num2str(numIterU),' res: ',num2str(errorU(end))])

                    numIterU = numIterU + 1;
                    uOld.fValues = obj.fem.uFun.fValues;
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
            obj.boundaryConditions = bcSet;
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


        %% %%%%%%%%%%%%%%%%%%%%%% ELASTIC EQUATION %%%%%%%%%%%%%%%%%%%%%%%% %%
        % function computeFEM(obj)
        %     quadOrder = 'QUADRATIC';
        %     quad = Quadrature.set(obj.mesh.type);
        %     quad.computeQuadrature(quadOrder);
        % 
        %     s.mesh = obj.mesh;
        %     s.type = 'ELASTIC';
        %     s.scale = 'MACRO';
        %     s.dim = '2D';
        %     s.material = obj.materialPhaseField.setMaterial(obj.phi,'Interpolated');
        %     s.boundaryConditions = obj.boundaryConditions;
        %     s.interpolationType = 'LINEAR';
        %     s.quadratureOrder = quadOrder;
        %     s.solverType = 'REDUCED';
        %     s.solverMode = 'DISP';
        %     obj.fem = PhysicalProblem.create(s);
        %     obj.fem.solve();
        % end

        function computeFEM(obj)
            LHS = obj.createInternalEnergyStiffnessMatrix();
            RHS = obj.computeElasticResidual();
            obj.deltaU = LHS\RHS;
        end

        function res = computeElasticResidual(obj)
            Fint = obj.createInternalEnergyElasticForceVector();
            Fext = obj.createExternalWorkForceVector();
            res = Fint + Fext;
        end

        function K = createInternalEnergyStiffnessMatrix(obj)
            [Huu,~] = obj.functional.energy.computeHessian(obj.u,obj.phi,'QUADRATIC');
            K = Huu;
        end

        function Fint = createInternalEnergyElasticForceVector(obj)
            [Ju,~] = obj.functional.energy.computeGradient(obj.u,obj.phi,'LINEAR');
            Fint = Ju;
        end

        function Fext = createExternalWorkForceVector(obj)
            fExt = obj.boundaryConditions.pointloadFun;
            [J] = obj.functional.extWork.computeGradient(obj.u,fExt,'LINEAR');
            Fext = J;
        end

        %% %%%%%%%%%%%%%%%% PHASE-FIELD EQUATION %%%%%%%%%%%%%%%%%%%%%%%% %%
        % Matrix equation
        function computePhaseField(obj)
            LHS = obj.computePhaseFieldLHS();
            RHS = obj.computePhaseFieldResidual();
            obj.deltaPhi = -LHS\RHS;
            %obj.deltaPhi = -obj.tau*RHS;
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
            K = obj.functional.nonLocalDamage.computeHessian(obj.phi,'CONSTANT');
        end

        %%% RHS %%%
        % Internal energy force vector
        function Fi = createInternalEnergyForceVector(obj)
            [~,Jphi] = obj.functional.energy.computeGradient(obj.fem.uFun,obj.phi,'LINEAR');
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
            fExt = obj.boundaryConditions.pointloadFun;
            totVal = obj.functional.extWork.computeFunction(obj.u,fExt,'LINEAR');
        end

        %% %%%%%%%%%%%%%%%%%% AUXILIARY METHODS %%%%%%%%%%%%%%% %%



        function totReact = computeReaction(obj)
            UpSide  = max(obj.mesh.coord(:,2));
            isInUp = abs(obj.mesh.coord(:,2)-UpSide)< 1e-12;
            nodes = 1:obj.mesh.nnodes;

            % bMesh = obj.mesh.createBoundaryMesh{4};
            % s.mesh = bMesh.mesh;
            % s.fValues = obj.fem.reactions(2*nodes(isInUp));
            % R = P1Function(s);
            % 
            % q.mesh = bMesh.mesh;
            % q.quadType = 'LINEAR';
            % q.type = 'Function';
            % int = Integrator.create(q);
            % totReact = int.compute(R);

            totReact = sum(obj.fem.reactions(2*nodes(isInUp)));

            % LeftSide  = min(obj.mesh.coord(:,1));
            % isInLeft = abs(obj.mesh.coord(:,1)-LeftSide)< 1e-12;
            % nodes = 1:obj.mesh.nnodes;
            % totReact = sum(obj.fem.reactions(2*nodes(isInLeft)));

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
            s.fValues = obj.fem.uFun.fValues - f.fValues;
            f = P1Function(s);
            int = Integrator.create('ScalarProduct',obj.mesh,'CONSTANT');
            error = sqrt(int.compute(f,f));

            error2 = norm(s.fValues);
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

            % obj.forceMat(step) = obj.fem.stressFun.fValues(2,1,1);  %% ONLY ONE ELEMENT
            obj.reactionMat(step) = obj.computeReaction();
            
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

            figure(100)
            plot(obj.displacementMat(1:step),obj.forceMat(1:step))
            title('Force-displacement diagram (Internal force)')
            grid on
            xlabel('Displacement [mm]')
            ylabel('Force [kN]')

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
