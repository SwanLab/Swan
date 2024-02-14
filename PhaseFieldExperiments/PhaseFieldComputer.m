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
        phaseField
        deltaPhi
        deltaU
    end

    properties (Access = private)
        energyFunc
        localDamageFunc
        nonLocalDamageFunc
        externalWorkFunc

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

            phaseFieldOld = obj.phaseField.fValues;
            Uold = P1Function.create(obj.mesh,2);

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
                        obj.phaseField.fValues = obj.phaseField.fValues + obj.deltaPhi;
                        obj.phaseField.fValues = min(max(phaseFieldOld, obj.phaseField.fValues),1);

                        obj.costFun(1,end+1) = obj.computeCostFunction();
                        obj.costFun(2,end) = 1;

                        errorPhi = abs(obj.costFun(1,end)-obj.costFun(1,end-1));
                        disp(['iterPhi: ',num2str(numIterP),' res: ',num2str(errorPhi)])

                        numIterP = numIterP + 1;
                    end
                    errorU(end+1) = obj.computeDisplacementError(Uold);
                    disp(['iterU: ',num2str(numIterU),' res: ',num2str(errorU(end))])

                    numIterU = numIterU + 1;
                    Uold.fValues = obj.fem.uFun.fValues;
                end
                phaseFieldOld = obj.phaseField.fValues;

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
            obj.phaseField               = cParams.initialPhaseField;
            obj.materialPhaseField       = cParams.materialPhaseField;
            obj.dissipationInterpolation = cParams.dissipationPhaseField;
            obj.l0 = cParams.l0;
            obj.constant = cParams.constant;

            obj.forceMat = zeros(1,obj.steps);
            obj.displacementMat = zeros(1,obj.steps);
            obj.damageMat = zeros(1,obj.steps);
            obj.stressStrainMat = zeros(2,obj.steps);
            obj.energyMat = zeros(4,obj.steps);
            obj.iterMat = zeros(2,obj.steps);
        end

        function createBoundaryConditions(obj,prescribedVal)
            s.mesh = obj.mesh;
            bc = BoundaryContionsForPhaseFieldCreator(s);
            bC = bc.create(prescribedVal);
            obj.boundaryConditions = bC;
        end

        function createFunctionals(obj)
            s.mesh = obj.mesh;
            s.materialPhaseField = obj.materialPhaseField;
            s.dissipationInterpolation = obj.dissipationInterpolation;
            s.constant = obj.constant;
            s.l0 = obj.l0;

            obj.energyFunc = ShFunc_InternalEnergy(s);
            obj.localDamageFunc = ShFunc_LocalDamage(s);
            obj.nonLocalDamageFunc = ShFunc_NonLocalDamage(s);
        end


        %% %%%%%%%%%%%%%%%%%%%%%% ELASTIC EQUATION %%%%%%%%%%%%%%%%%%%%%%%% %%
        function computeFEM(obj)
            quadOrder = 'QUADRATIC';
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature(quadOrder);
            obj.materialPhaseField.computeInterpolatedMaterial(obj.phaseField,quad);

            s.mesh = obj.mesh;
            s.type = 'ELASTIC';
            s.scale = 'MACRO';
            s.dim = '2D';
            s.material = obj.materialPhaseField.material;
            s.boundaryConditions = obj.boundaryConditions;
            s.interpolationType = 'LINEAR';
            s.quadratureOrder = quadOrder;
            s.solverType = 'REDUCED';
            s.solverMode = 'DISP';
            obj.fem = FEM.create(s);
            obj.fem.solve();
        end

        function computeFEM(obj)
            LHS = obj.createInternalEnergyStiffnessMatrix();
            RHS = obj.computeElasticResidual(obj);
            obj.deltaU = LHS\RHS;
        end


        function res = computeElasticResidual(obj)
            F = obj.createInternalEnergyElasticForceVector();
            res = F;
        end

        function K = createInternalEnergyStiffnessMatrix(obj)
            quadOrder = 'LINEAR';
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature(quadOrder);

            [Huu,~] = obj.energyFunc.computeHessian(obj.fem.uFun,obj.phaseField,quad);
            K = Huu;
        end

        function F = createInternalEnergyElasticForceVector(obj)
            quadOrder = 'LINEAR';
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature(quadOrder);
            [Ju,~] = obj.energyFunc.computeGradient(obj.fem.uFun,obj.phaseField,quad);
            F = Ju;
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
            quadOrder = 'QUADRATIC';
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature(quadOrder);

            [~,Hphiphi] = obj.energyFunc.computeHessian(obj.fem.uFun,obj.phaseField,quad);
            Mi = Hphiphi;
        end

        % Dissipation mass matrix
        function Md = createDissipationMassMatrix(obj)
            Md = obj.localDamageFunc.computeHessian(obj.phaseField,'LINEAR');
        end

        % Stiffness matrix
        function K = createStiffnessMatrix(obj)
            K = obj.nonLocalDamageFunc.computeHessian(obj.phaseField,'CONSTANT');
        end

        %%% RHS %%%
        % Internal energy force vector
        function Fi = createInternalEnergyForceVector(obj)
            quadOrder = 'LINEAR';
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature(quadOrder);
            [Ju,Jphi] = obj.energyFunc.computeGradient(obj.fem.uFun,obj.phaseField,quad);
            Fi = Jphi;
        end
        
        % Dissipation force vector
        function Fd = createDissipationForceVector(obj)
            Fd = obj.localDamageFunc.computeGradient(obj.phaseField,'LINEAR');     
        end
       
        % Force derivative vector
        function DF = createForceDerivativeVector(obj)
            quadOrder = 'LINEAR';
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature(quadOrder);
            DF = obj.nonLocalDamageFunc.computeGradient(obj.phaseField,quad);
        end

        %% %%%%%%%%%%%%%%%% COMPUTE TOTAL ENERGIES %%%%%%%%%%%%%%%%%%%%%%%% %%
        function totVal = computeTotalInternalEnergy(obj)
            quadOrder = 'QUADRATIC';
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature(quadOrder);
            totVal = obj.energyFunc.computeFunction(obj.fem.uFun,obj.phaseField,quad);
        end

        function totVal = computeTotalDissipationLocal(obj)
            totVal = obj.localDamageFunc.computeFunction(obj.phaseField,'QUADRATICMASS');
        end

        function totVal = computeTotalRegularizationTerm(obj)
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature('CONSTANT');
            totVal = obj.nonLocalDamageFunc.computeFunction(obj.phaseField,quad);
        end

        function totVal = computeTotalExternalWork(obj)
            % u = obj.fem.uFun;
            % forceValues = zeros(size(u.fValues));
            % pLoad = obj.boundaryConditions.pointload;
            % if isempty(pLoad)
            %     totVal = 0;
            % else 
            %     idx = sub2ind(size(forceValues),pLoad(:,1),pLoad(:,2));
            %     forceValues(idx) = pLoad(:,3);
            % 
            %     s.mesh = obj.mesh;
            %     s.fValues = forceValues;
            %     f = P1Function(s);
            % 
            %     q.mesh = obj.mesh;
            %     q.quadType = 'CONSTANT';
            %     q.type = 'ScalarProduct';
            %     int = Integrator.create(q);
            %     totVal = int.compute(u,f);

            pLoad = obj.boundaryConditions.pointloadFun;
            if isempty(pLoad)
                totVal = 0;
            else
                UpSide  = max(obj.mesh.coord(:,2));
                isInUp = abs(obj.mesh.coord(:,2)-UpSide)< 1e-12;
                nodes = 1:obj.mesh.nnodes;

                bMesh = obj.mesh.createBoundaryMesh{4};
                f = obj.boundaryConditions.pointloadFun;

                s.fValues = obj.fem.uFun.fValues(nodes(isInUp),2);
                u = P1Function(s);

                q.mesh = bMesh.mesh;
                q.quadType = 'LINEAR';
                q.type = 'ScalarProduct';
                int = Integrator.create(q);
                totVal = int.compute(u,f);
            end
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

            q.mesh = obj.mesh;
            q.quadType = 'CONSTANT';
            q.type = 'ScalarProduct';
            int = Integrator.create(q);
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

            obj.forceMat(step) = obj.fem.stressFun.fValues(2,1,1);  %% ONLY ONE ELEMENT
            obj.reactionMat(step) = obj.computeReaction();
            
            obj.displacementMat(step) = obj.bcVal(step);
            obj.damageMat(step) = max(obj.phaseField.fValues);
         
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

            obj.stressStrainMat(1,step);
            obj.stressStrainMat(2,step);

            obj.iterMat(1,step) = cParams.numIterU;              
            obj.iterMat(2,step) = cParams.numIterP;

            obj.computeEnergyIntegrals(step)
        end

        function printPlots(obj,step)
            figure(400)
            hold on
            obj.phaseField.plot;
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
            obj.phaseField.print(['PhaseFieldDamage_Step',step],'GiD');
        end

    end

end
