classdef PhaseFieldComputer < handle

    properties (Constant, Access = public)
        tolErrU = 1e-12;
        tolErrPhi = 1e-12;
        steps = 1000;
        tau = 1e2;
    end

    properties (Access = private)
        mesh
        boundaryConditions
        
        materialPhaseField
        dissipationInterpolation

        fem
        phaseField
        deltaPhi
    end

    properties (Access = private)
        Mi
        Md
        K
        Fi
        Fd
        DF
        Constant
        l0

        forceMat
        displacementMat
        damageMat
        stressStrainMat
        energyMat
        iterMat
    end

    methods (Access = public)

        function obj = PhaseFieldComputer(cParams)
            obj.init(cParams)

            obj.Constant = obj.materialPhaseField.Gc/(4*0.5);
            obj.l0 = 1*1e-1;
            
            phaseFieldOld = obj.phaseField.fValues;
            Uold = P1Function.create(obj.mesh,2);
            for i = 1:obj.steps
                disp([newline '%%%%%%%%%% STEP ',num2str(i),' %%%%%%%%%%%'])
                obj.createBoundaryConditions(i,obj.steps);
                numIterU = 1;
                errorU = 1;
                costFun = null(2,1);
                while (errorU(end) > obj.tolErrU) && (numIterU < 100)
                    obj.computeFEM();
                    costFun(1,end+1) = obj.computeCostFunction();
                    costFun(2,end) = 0;

                    numIterP = 1;
                    errorPhi = 1;
                    while (errorPhi > obj.tolErrPhi) && (numIterP < 100)
                        obj.solvePhaseFieldEquation();
                        obj.phaseField.fValues = obj.phaseField.fValues + obj.deltaPhi;
                        obj.phaseField.fValues = min(max(phaseFieldOld, obj.phaseField.fValues),1);

                        costFun(1,end+1) = obj.computeCostFunction();
                        costFun(2,end) = 1;

                        errorPhi = abs(costFun(1,end)-costFun(1,end-1));
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
                obj.printPlots(i);

                figure(500)
                for k = 2:size(costFun,2)
                    if costFun(2,k) == 0
                        points = [k   , costFun(1,k-1);
                                  k   , costFun(1,k)  ;
                                  k+1 , costFun(1,k)  ];
                        plot(points(:,1),points(:,2),'b')
                        hold on
                    elseif costFun(2,k) == 1
                        points = [k   , costFun(1,k-1);
                                  k   , costFun(1,k)  ;
                                  k+1 , costFun(1,k)  ];
                        plot(points(:,1),points(:,2),'r')
                        hold on
                    end

                end
                hold off
            end
        end

    end




    methods (Access = private)

        function init(obj,cParams)
            obj.mesh                     = cParams.mesh;
            obj.phaseField               = cParams.initialPhaseField;
            obj.materialPhaseField       = cParams.materialPhaseField;
            obj.dissipationInterpolation = cParams.dissipationPhaseField;

            obj.forceMat = zeros(1,obj.steps);
            obj.displacementMat = zeros(1,obj.steps);
            obj.damageMat = zeros(1,obj.steps);
            obj.stressStrainMat = zeros(2,obj.steps);
            obj.energyMat = zeros(4,obj.steps);
            obj.iterMat = zeros(2,obj.steps);
        end

        function createBoundaryConditions(obj,i,nIter)
            s.mesh = obj.mesh;
            bc = BoundaryContionsForPhaseFieldCreator(s);
            bC = bc.create(i,nIter);
            obj.boundaryConditions = bC;
        end


        %% %%%%%%%%%%%%%%%%%%%%%% ELASTIC EQUATION %%%%%%%%%%%%%%%%%%%%%%%% %%
        function computeFEM(obj)
            quadOrder = 'QUADRATIC';
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature(quadOrder);

            sMat.quadrature = quad;
            sMat.phi = obj.phaseField;
            sMat.derivative = 0;
            obj.materialPhaseField.computeMatInt(sMat);

            s.mesh = obj.mesh;
            s.type = 'ELASTIC';
            s.scale = 'MACRO';
            s.dim = '2D';
            s.material = obj.materialPhaseField.material;
            s.bc = obj.boundaryConditions;
            s.interpolationType = 'LINEAR';
            s.quadratureOrder = quadOrder;
            obj.fem = FEM.create(s);
            obj.fem.solve();
        end

        %% %%%%%%%%%%%%%%%% PHASE-FIELD EQUATION (LHS) %%%%%%%%%%%%%%%%%%%%%%%% %%
        % Internal energy mass matrix
        function createInternalEnergyMassMatrix(obj)
            DDenergyFun =  obj.createAbstractDerivativeEnergyFunction('LINEAR',2); 

            s.function = DDenergyFun;
            s.trial = P1Function.create(obj.mesh,1);
            s.test = P1Function.create(obj.mesh,1);
            s.mesh = obj.mesh;
            s.type = 'MassMatrixWithFunction';
            s.quadratureOrder = 'LINEAR';
            LHS = LHSintegrator.create(s);
            obj.Mi = LHS.compute(); 
        end

        % Dissipation mass matrix
        function createDissipationMassMatrix(obj)
            ddAlphaFun =  obj.createSecondDerivativeDissipationFunction();

            s.trial = P1Function.create(obj.mesh,1);
            s.test = P1Function.create(obj.mesh,1);
            s.function = ddAlphaFun;
            s.mesh = obj.mesh;
            s.type = 'MassMatrixWithFunction';
            s.quadratureOrder = 'LINEAR';
            LHS = LHSintegrator.create(s);
            obj.Md = LHS.compute(); 
        end

        function ddAlpha = createSecondDerivativeDissipationFunction(obj)
            s.mesh = obj.mesh;
            s.handleFunction = obj.dissipationInterpolation.ddfun;
            s.l2function = obj.phaseField;
            ddAlpha = CompositionFunction(s);
        end

        % Stiffness matrix
        function createStiffnessMatrix(obj)
            s.trial = P1Function.create(obj.mesh,1);
            s.test = P1Function.create(obj.mesh,1);
            s.mesh = obj.mesh;
            s.type = 'StiffnessMatrix';
            LHS = LHSintegrator.create(s);
            obj.K = LHS.compute();  
        end

        %% %%%%%%%%%%%%%%%% PHASE-FIELD EQUATION (RHS) %%%%%%%%%%%%%%%%%%%%%%%% %%
        % Internal energy force vector
        function createInternalEnergyForceVector(obj)
            DenergyFun =  obj.createAbstractDerivativeEnergyFunction('LINEAR',1); 
            test = P1Function.create(obj.mesh,1);
            
            s.mesh = obj.mesh;
            s.type = 'ShapeFunction';
            s.quadType = 'LINEAR';
            RHS = RHSintegrator.create(s);
            obj.Fi = RHS.compute(DenergyFun,test); 
        end
        
        % Dissipation force vector
        function createDissipationForceVector(obj)
            dAlphaFun =  obj.createFirstDerivativeDissipationFunction(); 
            test = P1Function.create(obj.mesh,1);

            s.mesh = obj.mesh;
            s.type = 'ShapeFunction';
            s.quadType = 'LINEAR';
            RHS = RHSintegrator.create(s);
            obj.Fd = RHS.compute(dAlphaFun, test);     
        end

        function dAlpha = createFirstDerivativeDissipationFunction(obj)
            s.mesh = obj.mesh;
            s.handleFunction = obj.dissipationInterpolation.dfun;
            s.l2function = obj.phaseField;
            dAlpha = CompositionFunction(s);
        end  
       
        % Force derivative vector
        function createForceDerivativeVector(obj)
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature('LINEAR');

            PhiGradient = obj.phaseField.computeGradient(quad);
            test = P1Function.create(obj.mesh,1);

            s.quadratureOrder = 'LINEAR';
            s.mesh = obj.mesh;
            s.type = 'ShapeDerivative';
            RHS = RHSintegrator.create(s);
            forceVector = RHS.compute(PhiGradient, test); 
            obj.DF = forceVector.fValues;
        end

        % Matrix equation
        function solvePhaseFieldEquation(obj)
            obj.createInternalEnergyMassMatrix();
            obj.createDissipationMassMatrix();
            obj.createStiffnessMatrix();
           

            LHS = obj.Mi + (obj.Constant/obj.l0)*obj.Md + (obj.Constant*obj.l0)*obj.K;
            RHS = obj.computeResidual();
            obj.deltaPhi = -LHS\RHS;
            %obj.deltaPhi = -obj.tau*RHS;
        end

        function res = computeResidual(obj)
            obj.createInternalEnergyForceVector();
            obj.createDissipationForceVector();
            obj.createForceDerivativeVector();            
            res = (obj.Fi + (obj.Constant/obj.l0)*obj.Fd + (obj.Constant*obj.l0)*obj.DF);
        end

        %% %%%%%%%%%%%%%%%% COMPUTE TOTAL ENERGIES %%%%%%%%%%%%%%%%%%%%%%%% %%
        function totVal = computeTotalEnergy(obj)
            e = obj.createAbstractDerivativeEnergyFunction('QUADRATICMASS',0);

            q.mesh = obj.mesh;
            q.quadType = 'QUADRATICMASS';
            q.type = 'Function';
            int = Integrator.create(q);
            totVal = int.compute(e);
        end

        function totVal = computeTotalDissipationLocal(obj)
            alphaFun = obj.createDissipationFunction();

            q.mesh = obj.mesh;
            q.quadType = 'QUADRATICMASS';
            q.type = 'Function';
            int = Integrator.create(q);
            totVal = (obj.Constant/obj.l0)*int.compute(alphaFun);
        end

        function dAlpha = createDissipationFunction(obj)
            s.mesh = obj.mesh;
            s.handleFunction = obj.dissipationInterpolation.dfun;
            s.l2function = obj.phaseField;
            dAlpha = CompositionFunction(s);
        end  

        function totVal = computeTotalRegularizationTerm(obj)
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature('CONSTANT');

            PhiGradient = obj.phaseField.computeGradient(quad);
            GradGrad = sum(PhiGradient.fValues.^2);
            
            s.fValues = GradGrad;
            s.quadrature = quad;
            s.mesh = obj.mesh;
            GradGradFun = FGaussDiscontinuousFunction(s);
            
            q.mesh = obj.mesh;
            q.quadType = 'CONSTANT';
            q.type = 'Function';
            int = Integrator.create(q);
            totVal = 0.5*(obj.Constant*obj.l0)*int.compute(GradGradFun);
        end

        function totVal = computeTotalExternalWork(obj)
            u = obj.fem.uFun;
            forceValues = zeros(size(u.fValues));
            pLoad = obj.boundaryConditions.pointload;
            if isempty(pLoad)
                totVal = 0;
            else 
                idx = sub2ind(size(forceValues),pLoad(:,1),pLoad(:,2));
                forceValues(idx) = pLoad(:,3);

                s.mesh = obj.mesh;
                s.fValues = forceValues;
                f = P1Function(s);

                q.mesh = obj.mesh;
                q.quadType = 'CONSTANT';
                q.type = 'ScalarProduct';
                int = Integrator.create(q);
                totVal = int.compute(u,f);
            end

        end
        %% %%%%%%%%%%%%%%%%%% AUXILIARY METHODS %%%%%%%%%%%%%%% %%
        function energyVal = computeEnergyField(obj,e,C)
            nstre = size(e.fValues,1);
            nGauss = size(e.fValues,2);
            nelem = size(e.fValues,3);
            energyVal = zeros(1,nGauss,nelem);
            for iStre = 1:nstre
                for jStre=1:nstre
                    for iGauss=1:nGauss
                        eI = squeeze(e.fValues(iStre,iGauss,:));
                        eJ = squeeze(e.fValues(jStre,iGauss,:));
                        Cij = squeeze(C(iStre,jStre,:,iGauss));
                        eStre(1,1,:) = (eI.*Cij.*eJ)';
                        energyVal(1,iGauss,:) = energyVal(1,iGauss,:) + eStre;
                    end
                end
            end
            energyVal = 0.5*energyVal;
        end

        function energyFun = createAbstractDerivativeEnergyFunction(obj,quadOrder,deriv)
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature(quadOrder);

            e = obj.fem.uFun.computeSymmetricGradient(quad);
            e.applyVoigtNotation();

            s.quadrature = quad;
            s.phi = obj.phaseField;
            s.derivative = deriv;
            obj.materialPhaseField.computeMatInt(s);
            C = obj.materialPhaseField.material.C;
            DDenergyVal = obj.computeEnergyField(e,C);

            s.fValues = DDenergyVal;
            s.quadrature = quad;
            s.mesh = obj.mesh;
            energyFun = FGaussDiscontinuousFunction(s);
        end
        

        function totVal = computeIntTotalForce(obj)
            stress = obj.fem.stressFun;

            q.mesh = obj.mesh;
            q.quadType = 'QUADRATIC';
            q.type = 'Function';
            int = Integrator.create(q);
            Vol = obj.mesh.computeVolume();

            totVal = int.compute(stress)/Vol;
        end

        function e = computeCostFunction(obj)
            eW = obj.computeTotalExternalWork();
            eT = obj.computeTotalEnergy();
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

        %% %%%%%%%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%%% %%

        function saveData(obj,cParams)
            step = cParams.step;

            obj.forceMat(step) = obj.computeIntTotalForce();
            obj.displacementMat(step) = max(abs(obj.fem.uFun.fValues(:,2)));
            obj.damageMat(step) = max(obj.phaseField.fValues);

            obj.energyMat(1,step) = obj.computeTotalExternalWork();
            obj.energyMat(2,step) = obj.computeTotalEnergy();
            obj.energyMat(3,step) = obj.computeTotalDissipationLocal();
            obj.energyMat(4,step) = obj.computeTotalRegularizationTerm();

            obj.stressStrainMat(1,step);
            obj.stressStrainMat(2,step);

            obj.iterMat(1,step) = cParams.numIterU;              
            obj.iterMat(2,step) = cParams.numIterP;
        end

        function printPlots(obj,step)
            figure(100)
            plot(obj.displacementMat(1:step),obj.forceMat(1:step))
            title('Force-displacement diagram')
            xlabel('Displacement [mm]')
            ylabel('Force [kN]')

            figure(101)
            plot(obj.displacementMat(1:step),obj.damageMat(1:step))
            title('Damage-displacement diagram')
            xlabel('Displacement [mm]')
            ylabel('Damage [-]')

            figure(200)
            hold on
            plot(obj.energyMat(1,1:step),'Color',"#0072BD")
            plot(obj.energyMat(2,1:step),'Color',"#D95319")
            plot(obj.energyMat(3,1:step),'Color',"#EDB120")
            plot(obj.energyMat(4,1:step),'Color',"#7E2F8E")
            title('Energy values at each step')
            legend('External Work', ...
                   'Internal Energy', ...
                   'Local surface energy', ...
                   'Non-local surface energy')
            xlabel('Steps [-]')
            ylabel('Energy [J]')
            
            figure(300)
            hold on
            plot(obj.iterMat(1,1:step),'b')
            plot(obj.iterMat(2,1:step),'r')
            title('Iterations needed')
            legend('U','phi')
            xlabel('Step')
            ylabel('Iterations')
            
            figure(400)
            obj.phaseField.plot;
            colorbar
            clim([0 1])
            drawnow
            title(['Damage at step ',num2str(step)])
        end

        function printResults(obj,step)
            obj.fem.print(['PhaseFieldFEM_Step',step],'GiD');
            obj.phaseField.print(['PhaseFieldDamage_Step',step],'GiD');
        end

    end

end
