classdef PhaseFieldComputer < handle

    properties (Constant, Access = public)
        tolErrU = 1e-12;
        tolErrPhi = 1e-12;
        Gc = 5e-3;
        fc = 1;

        E = 210;
        nu = 0.3;
    end

    properties (Access = private)
        mesh
        quadrature
        boundaryConditions
        materialInterpolation
        dissipationInterpolation
        phaseField
        deltaPhi
        fem
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
            obj.createQuadrature();
            obj.createPhaseField();
            obj.createMaterialInterpolation();
            obj.createDissipationInterpolation();

            %obj.l0 = (27/256)*(1*obj.Gc/obj.fc^2);
            obj.Constant = obj.Gc/(4*0.5);

            obj.l0 = 10*1e-1;

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
                disp(['Setp Number: ',num2str(i)])
                while (errorU > obj.tolErrU) && (numIterU < 100)
                    obj.computeFEM();
        %            obj.fem.uFun.plot()
                    errorPhi = 1;
                    numIterP = 1;
                    while (errorPhi > obj.tolErrPhi) && (numIterP < 100)

                        obj.solvePhaseFieldEquation();
                        %obj.computeTotalEnergy()

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
                plot(ForceDisplacement(2,:),ForceDisplacement(1,:))
                figure(101)
                plot(ForceDisplacement(2,:),PhaseField)
                %obj.fem.uFun.plot;
                
                 
                 obj.phaseField.plot;                
                 colorbar
                 caxis([0 1]);
                 %clim([0 1])
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
            obj.mesh               = cParams.mesh; 
        end

        function createQuadrature(obj)
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature('LINEAR');
            obj.quadrature = quad;
        end

        function createBoundaryConditions(obj,i,nIter)
            s.mesh = obj.mesh;
            bc = BoundaryContionsForPhaseFieldCreator(s);
            bC = bc.create(i,nIter);
            obj.boundaryConditions = bC;
        end

        function createPhaseField(obj)
            % sLS.type       = 'rectangleInclusion';
            % sLS.mesh       = obj.mesh;
            % sLS.ndim       = 2;
            % sLS.widthH = 1;
            % sLS.widthV = 2;
            % sLS.coord      = obj.mesh.coord;
            % ls = LevelSetCreator.create(sLS);
            % dmg = ls.getValue();
            % obj.phi = dmg;
            % scatter3(obj.mesh.coord(:,1),obj.mesh.coord(:,2),obj.phi);

            xmax = max(obj.mesh.coord(:,1));
            sAF.fHandle = @(x) x(1,:,:)-x(1,:,:);
            %sAF.fHandle = @(x) x(1,:,:)/xmax;
            %sAF.fHandle = @(x) (x(1,:,:).^2)/xmax;
            sAF.ndimf   = 1;
            sAF.mesh    = obj.mesh;
            xFun = AnalyticalFunction(sAF);

            phi = xFun.project('P1');
            %phi.fValues = 0.
            obj.phaseField = phi;
        end

        function createMaterialInterpolation(obj)
            c.typeOfMaterial = 'ISOTROPIC';
            c.interpolation = 'PhaseFieldI';
            c.nElem = obj.mesh.nelem;
            c.dim = '2D';
            c.constitutiveProperties.rho_plus = 1;
            c.constitutiveProperties.rho_minus = 0;
            c.constitutiveProperties.E_plus = obj.E;
            c.constitutiveProperties.E_minus = 1e-3;
            c.constitutiveProperties.nu_plus = obj.nu;
            c.constitutiveProperties.nu_minus = 1/3;

            matInt = MaterialInterpolation.create(c);
            obj.materialInterpolation = matInt;
        end

        function createDissipationInterpolation(obj)
            c.typeOfMaterial = 'ISOTROPIC';
            c.interpolation = 'PhaseFieldD';

            disInt = MaterialInterpolation.create(c);
            obj.dissipationInterpolation = disInt;
        end


        %% %%%%%%%%%%%%%%%%%%%%%% ELASTIC EQUATION %%%%%%%%%%%%%%%%%%%%%%%% %%
        function mat = createMaterial(obj)
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature('QUADRATIC');

            phiV = obj.phaseField.evaluate(quad.posgp);
            phiV = permute(phiV,[1 3 2]); 
            matInt  = obj.materialInterpolation.computeMatProp(squeezeParticular(phiV,1));
            s.ptype = 'ELASTIC';
            s.pdim  = '2D';
            s.nelem = obj.mesh.nelem;
            s.mesh  = obj.mesh;
            s.kappa = matInt.kappa;
            s.mu    = matInt.mu;
            mat = Material.create(s);
            mat.compute(s);
        end

        function computeFEM(obj)
            s.mesh = obj.mesh;
            s.type = 'ELASTIC';
            s.scale = 'MACRO';
            s.material = obj.createMaterial();
            s.dim = '2D';
            s.bc = obj.boundaryConditions;
            s.interpolationType = 'LINEAR';
            s.quadratureOrder = 'QUADRATIC';
            obj.fem = FEM.create(s);
            obj.fem.solve();
        end

        %% %%%%%%%%%%%%%%%% PHASE-FIELD EQUATION (LHS) %%%%%%%%%%%%%%%%%%%%%%%% %%

        % Internal energy mass matrix
        function DDenergy = createFGaussDDEnergyFunction(obj)
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature('LINEAR');

            e = obj.fem.uFun.computeSymmetricGradient(quad);
            e.applyVoigtNotation();
            phiV = obj.phaseField.evaluate(quad.posgp);
            phiV = permute(phiV,[1 3 2]); 
            matInt  = obj.materialInterpolation.computeDDMatProp(squeezeParticular(phiV,1));

            DDenergyV = obj.computeEnergyField(e,matInt.ddmu,matInt.ddkappa);

            s.fValues = DDenergyV;
            s.quadrature = quad;
            s.mesh = obj.mesh;
            DDenergy = FGaussDiscontinuousFunction(s);

        end

     
        function createEnergyMassMatrix(obj)
            DDenergyFun =  obj.createFGaussDDEnergyFunction(); 
            s.trial = P1Function.create(obj.mesh,1);
            s.test = P1Function.create(obj.mesh,1);
            s.function = DDenergyFun;
            s.mesh = obj.mesh;
            s.type = 'MassMatrixWithFunction';
            s.quadratureOrder = 'LINEAR';
            LHS = LHSintegrator.create(s);
            obj.Mi = LHS.compute(); 
        end
        
        % Dissipation mass matrix
        function DDalpha = createFGaussDDDissipationFunction(obj)
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature('LINEAR');

            phiV     = obj.phaseField.evaluate(quad.posgp);
            DDalphaV = obj.dissipationInterpolation.computeDDAlphaProp(phiV);
            s.fValues = DDalphaV;
            s.quadrature = obj.quadrature;
            s.mesh = obj.mesh;
            DDalpha = FGaussDiscontinuousFunction(s);
        end

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

        % Stiffness matrix
        function createStiffnessMatrix(obj)
            s.trial = P1Function.create(obj.mesh,1);
            s.test = P1Function.create(obj.mesh,1);
            s.mesh = obj.mesh;
            s.type = 'StiffnessMatrix';
            LHS = LHSintegrator.create(s);
            obj.K = LHS.compute();  
        end

        %% PHASE-FIELD EQUATION (LHS)
        % Internal energy force vector
        function Denergy = createFGaussDEnergyFunction(obj)
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature('LINEAR');

            e = obj.fem.uFun.computeSymmetricGradient(quad);
            e.applyVoigtNotation();
            phiV = obj.phaseField.evaluate(quad.posgp);
            phiV = permute(phiV,[1 3 2]); 
            matInt  = obj.materialInterpolation.computeDMatProp(squeezeParticular(phiV,1));

            DenergyV = obj.computeEnergyField(e,matInt.dmu,matInt.dkappa);

            s.fValues = DenergyV;
            s.quadrature = quad;
            s.mesh = obj.mesh;
            Denergy = FGaussDiscontinuousFunction(s);
            %Denergy.plot();
            %title('DEnergy')
        end

        function energyVal = computeEnergyField(obj,e,mu,kappa)

            s.ptype = 'ELASTIC';
            s.pdim  = '2D';
            s.nelem = obj.mesh.nelem;
            s.mesh  = obj.mesh;
            s.kappa = kappa;
            s.mu    = mu;
            mat = Material.create(s);
            mat.compute(s);

            nstre = size(e.fValues,1);
            nGauss = size(e.fValues,2);
            nelem = size(e.fValues,3);
            energyVal = zeros(1,nGauss,nelem);
            for iStre = 1:nstre
                for jStre=1:nstre
                    for iGauss=1:nGauss
                        eI = squeeze(e.fValues(iStre,iGauss,:));
                        eJ = squeeze(e.fValues(jStre,iGauss,:));
                        ddCij = squeeze(mat.C(iStre,jStre,:,iGauss));
                        eStre(1,1,:) = (eI.*ddCij.*eJ)';
                        energyVal(1,iGauss,:) = energyVal(1,iGauss,:) + eStre;
                    end
                end
            end
            energyVal = 0.5*energyVal;
        end

        function createEnergyForceVector(obj)
            DenergyFun =  obj.createFGaussDEnergyFunction(); 
            test = P1Function.create(obj.mesh,1);
            s.mesh = obj.mesh;
            s.type = 'ShapeFunction';

            RHS = RHSintegrator.create(s);
            obj.Fi = RHS.compute(DenergyFun,test); 
        end
        
        % Dissipation force vector
        function Dalpha = createFGaussDDissipationFunction(obj)
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature('LINEAR');

            phiV     = obj.phaseField.evaluate(quad.posgp);
            DalphaV = obj.dissipationInterpolation.computeDAlphaProp(phiV);
            s.fValues = DalphaV;
            s.quadrature = obj.quadrature;
            s.mesh = obj.mesh;
            Dalpha = FGaussDiscontinuousFunction(s);
        end   

        function createDissipationForceVector(obj)
            DalphaFun =  obj.createFGaussDDissipationFunction(); 
            test = P1Function.create(obj.mesh,1);
            s.mesh = obj.mesh;
            s.type = 'ShapeFunction';
            s.quadratureOrder = obj.quadrature.order;
            RHS = RHSintegrator.create(s);
            obj.Fd = RHS.compute(DalphaFun, test);     
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

        % Matrix euqation
        function solvePhaseFieldEquation(obj)
            obj.createEnergyMassMatrix();
            obj.createDissipationMassMatrix();
            obj.createStiffnessMatrix();
            obj.createEnergyForceVector();
            obj.createDissipationForceVector();
            obj.createForceDerivativeVector();
            LHS = obj.Mi + (obj.Constant/obj.l0)*obj.Md + (obj.Constant*obj.l0)*obj.K;
            RHS = obj.computeResidual();
            obj.deltaPhi = LHS\RHS;
        end

        function res = computeResidual(obj)
            res = -(obj.Fi + (obj.Constant/obj.l0)*obj.Fd + (obj.Constant*obj.l0)*obj.DF);
        end

        function totVal = computeTotalEnergy(obj)
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature('QUADRATIC');

            phiV = obj.phaseField.evaluate(quad.posgp);
            phiV = permute(phiV,[1 3 2]); 
            matInt  = obj.materialInterpolation.computeMatProp(squeezeParticular(phiV,1));
            s.ptype = 'ELASTIC';
            s.pdim  = '2D';
            s.nelem = obj.mesh.nelem;
            s.mesh  = obj.mesh;
            s.kappa = matInt.kappa;
            s.mu    = matInt.mu;
            mat = Material.create(s);
            mat.compute(s);
            e = obj.fem.uFun.computeSymmetricGradient(quad);
            e.applyVoigtNotation();
            energyV = obj.computeEnergyField(e,matInt.mu,matInt.kappa);

            s.fValues = energyV;
            s.quadrature = quad;
            s.mesh = obj.mesh;
            energyFun = FGaussDiscontinuousFunction(s);
            %energyFun.plot()
            %title('Energy')

            q.mesh = obj.mesh;
            q.quadType = 'QUADRATIC';
            q.type = 'InternalEnergy';
            int = Integrator.create(q);
            totVal = int.compute(e,mat.C);
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
