classdef PhaseFieldComputer < handle

    properties (Constant, Access = public)
        tolErrU = 1e-12;
        tolErrPhi = 1e-12;
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
    end

    methods (Access = public)

        function obj = PhaseFieldComputer(cParams)
            obj.init(cParams)

            %obj.l0 = (27/256)*(1*obj.material.Gc/obj.material.fc^2);
            obj.Constant = obj.materialPhaseField.Gc/(4*0.5);

            obj.l0 = 1*1e-1;

            niter = 1000;
            Energy = zeros(4,niter);
            ForceDisplacement = zeros(2,niter);
            Iterations = zeros(2,niter);
            F = zeros(3,niter);
            PhaseField = zeros(1,niter);
            for i = 1:niter
                
                obj.createBoundaryConditions(i,niter);
                errorU = 1;
                Uold = zeros(obj.mesh.nnodes,obj.mesh.ndim);
                numIterU = 1;
                numIterP = 1;

                disp([newline '%%%%%%%%%% NEW STEP %%%%%%%%%%%'])
                disp(['Step Number: ',num2str(i)])
                while (errorU > obj.tolErrU) && (numIterU < 100)
                    obj.computeFEM();
                    %obj.fem.uFun.plot()
                    errorPhi = 1;
                    numIterP = 1;
                    while (errorPhi > obj.tolErrPhi) && (numIterP < 100)
                        obj.solvePhaseFieldEquation();

                        obj.phaseField.fValues = obj.phaseField.fValues + obj.deltaPhi;
                        res = obj.computeResidual();
                        errorPhi = norm(res);
                        disp(['iterPhi: ',num2str(numIterP),' res: ',num2str(errorPhi)])
                        numIterP = numIterP + 1;
                        %obj.phaseField.plot;  
                    end
                    errorU = norm(obj.fem.uFun.fValues - Uold);
                    disp(['iterU: ',num2str(numIterU),' res: ',num2str(errorU)])
                    Uold = obj.fem.uFun.fValues;
                    numIterU = numIterU + 1;
                end

                F(1,i) = mean(obj.Fi);
                F(2,i) = mean(obj.Fd);
                F(3,i) = mean(obj.DF);

                Energy(1,i) = obj.computeTotalExternalWork();
                Energy(2,i) = obj.computeTotalEnergy();
                Energy(3,i) = obj.computeTotalDissipationLocal();
                Energy(4,i) = obj.computeTotalDissipationNonLocal();
                Iterations(1,i) = numIterU;
                Iterations(2,i) = numIterP;
                ForceDisplacement(1,i) = obj.computeIntTotalForce();
                ForceDisplacement(2,i) = max(abs(obj.fem.uFun.fValues(:,2)));
                PhaseField(1,i) = max(obj.phaseField.fValues);


                figure(100)
                plot(ForceDisplacement(2,1:i),ForceDisplacement(1,1:i))
                figure(101)
                plot(ForceDisplacement(2,1:i),PhaseField(1:i))
                %obj.fem.uFun.plot;
                
                 
                 obj.phaseField.plot;                
                 colorbar
                 clim([0 1])
                 drawnow
                % title('Final phase field')
                % %obj.fem.print(['Example',num2str(i)],'GiD')
                %obj.phaseField.print(['Example',num2str(i),'Phi'],'GiD')
            end
            figure
            plot(F(1,:))
            hold on
            plot(F(2,:))
            hold on
            plot(F(3,:))
            legend('Fi','Fd','DF')


            figure
            plot(Energy(1,:))
            hold on
            plot(Energy(2,:))
            hold on
            plot(Energy(3,:))
            hold on
            plot(Energy(4,:))
            legend('External Work','Internal Energy','Local surface energy','Non-local surface energy')

            figure
            plot(ForceDisplacement(2,:),ForceDisplacement(1,:))
            figure
            plot(ForceDisplacement(2,:))

            figure
            plot(Iterations(1,:))
            hold on
            plot(Iterations(2,:))
            legend('U','phi')
        end

    end




    methods (Access = private)

        function init(obj,cParams)
            obj.mesh                     = cParams.mesh;
            obj.phaseField               = cParams.initialPhaseField;
            obj.materialPhaseField       = cParams.materialPhaseField;
            obj.dissipationInterpolation = cParams.dissipationPhaseField;
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
            DDenergyFun =  obj.createEnergyFunction('LINEAR',2); 

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
            DDalphaFun =  obj.createFGaussDDDissipationFunction(); 
            s.trial = P1Function.create(obj.mesh,1);
            s.test = P1Function.create(obj.mesh,1);
            s.function = DDalphaFun;
            s.mesh = obj.mesh;
            s.type = 'MassMatrixWithFunction';
            s.quadratureOrder = 'LINEAR';
            LHS = LHSintegrator.create(s);
            obj.Md = LHS.compute(); 
        end

        function DDalpha = createFGaussDDDissipationFunction(obj)
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature('LINEAR');

            phiV     = obj.phaseField.evaluate(quad.posgp);
            DDalphaV = obj.dissipationInterpolation.computeDDAlphaProp(phiV);
            s.fValues = DDalphaV;
            s.quadrature = quad;
            s.mesh = obj.mesh;
            DDalpha = FGaussDiscontinuousFunction(s);
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
            DenergyFun =  obj.createEnergyFunction('LINEAR',1); 
            test = P1Function.create(obj.mesh,1);
            
            s.mesh = obj.mesh;
            s.type = 'ShapeFunction';
            RHS = RHSintegrator.create(s);
            obj.Fi = RHS.compute(DenergyFun,test); 
        end
        
        % Dissipation force vector
        function createDissipationForceVector(obj)
            DalphaFun =  obj.createFGaussDDissipationFunction(); 
            test = P1Function.create(obj.mesh,1);

            s.mesh = obj.mesh;
            s.type = 'ShapeFunction';
            s.quadratureOrder = 'LINEAR';
            RHS = RHSintegrator.create(s);
            obj.Fd = RHS.compute(DalphaFun, test);     
        end

        function Dalpha = createFGaussDDissipationFunction(obj)
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature('LINEAR');

            phiV     = obj.phaseField.evaluate(quad.posgp);
            DalphaV = obj.dissipationInterpolation.computeDAlphaProp(phiV);
            s.fValues = DalphaV;
            s.quadrature = quad;
            s.mesh = obj.mesh;
            Dalpha = FGaussDiscontinuousFunction(s);
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
            obj.createInternalEnergyForceVector();
            obj.createDissipationForceVector();
            obj.createForceDerivativeVector();

            LHS = obj.Mi + (obj.Constant/obj.l0)*obj.Md + (obj.Constant*obj.l0)*obj.K;
            RHS = obj.computeResidual();
            obj.deltaPhi = LHS\RHS;
        end

        function res = computeResidual(obj)
            res = -(obj.Fi + (obj.Constant/obj.l0)*obj.Fd + (obj.Constant*obj.l0)*obj.DF);
        end
        %% %%%%%%%%%%%%%%%% COMPUTE TOTAL ENERGIES %%%%%%%%%%%%%%%%%%%%%%%% %%
        function totVal = computeTotalEnergy(obj)
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature('QUADRATIC');

            e = obj.fem.uFun.computeSymmetricGradient(quad);
            e.applyVoigtNotation();

            obj.materialPhaseField.computeMatIso(quad);
            C = obj.materialPhaseField.material.C;

            q.mesh = obj.mesh;
            q.quadType = 'QUADRATIC';
            q.type = 'InternalEnergy';
            int = Integrator.create(q);
            totVal = int.compute(e,C);
        end

        function totVal = computeTotalDissipationLocal(obj)
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature('QUADRATIC');

            phiV = obj.phaseField.evaluate(quad.posgp);
            aValues = obj.dissipationInterpolation.computeAlphaProp(phiV);
           
            s.fValues = aValues;
            s.quadrature = quad;
            s.mesh = obj.mesh;
            alpha = FGaussDiscontinuousFunction(s);

            q.mesh = obj.mesh;
            q.quadType = 'QUADRATIC';
            q.type = 'Function';
            int = Integrator.create(q);
            totVal = (obj.Constant/obj.l0)*int.compute(alpha);
        end

        function totVal = computeTotalDissipationNonLocal(obj)
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

        function energyFun = createEnergyFunction(obj,quadOrder,deriv)
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
    end

end
