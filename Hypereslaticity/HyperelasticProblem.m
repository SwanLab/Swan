classdef HyperelasticProblem < handle

    properties (Access = public)
        uFun
        rFun
    end

    properties (Access = private)
        print
        bcCase
        fileNameOut
        meshType
    end

    properties (Access = private)
        mesh
        matProp, material
        neohookeanFun, linearElasticityFun
        monitor
        boundaryConditions

        bcApplier
        freeDofs, constrainedDofs
    end

    methods (Access = public)

        function obj = HyperelasticProblem(cParams)
            obj.init(cParams);
            obj.createMesh();
            obj.createBoundaryConditions();
            obj.createMaterial();
            obj.createDisplacementFun();
            obj.createFunctionals();
            obj.setMonitoring();
        end

        function solve(obj)
            u  = obj.uFun;
            iter = 1; nSteps = length(obj.boundaryConditions.bcValues);
            for iStep = 1:nSteps
                [u,bc] = obj.updateBoundaryConditions(u);

                %intEnergy(iter) = obj.neohookeanFun.compute(u);
                intEnergy(iter) = obj.linearElasticityFun.compute(u);
                Res = obj.computeResidual(u);
                resi(iter) = norm(Res);

                hasNotConverged = 1;
                while hasNotConverged
                    uVal = obj.reshapeToVector(u.fValues);
                    [hess, KR] = obj.computeHessian(u);
                    Res  = obj.computeResidual(u);
                    incU = hess\(-Res);
                    uVal(obj.freeDofs) = uVal(obj.freeDofs) + incU;
                    u.setFValues(obj.reshapeToMatrix(uVal));
                    reacFun = obj.computeReactions(KR,u);

                    iter = iter+1;
                    intEnergy(iter) = obj.neohookeanFun.compute(u);
                    Res        = obj.computeResidual(u);
                    resi(iter) = norm(Res);

                    hasNotConverged = obj.hasNotConverged(intEnergy);
                end
                normDisp(iStep) = Norm(u,'L2');
                normReac(iStep) = Norm(reacFun,'L2');
                obj.printFile(iStep, u, reacFun);
                nIterPerStep(iStep) = obj.computeNumberOfIterations(iter,nIterPerStep,iStep);
                obj.plotStep(intEnergy,nIterPerStep,iStep, normDisp,normReac);
                obj.rFun = reacFun;
            end
        end
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.print = cParams.printing;
            obj.bcCase  = cParams.bcCase;
            obj.fileNameOut = cParams.fileName;
            obj.meshType  = cParams.meshGen;
        end

        function createMesh(obj)
            switch obj.meshType
                case {'Hole', 'HoleDirich'}
                    IM = Mesh.createFromGiD('holeMeshQuad.m');
                    obj.mesh = IM;
                case {'Bending', 'Traction'}
                    obj.mesh = UnitQuadMesh(20,20);
                case {'Metamaterial'}
                    load('NegativePoissonMesh.mat','NegPoissMesh');
                    s.coord  = NegPoissMesh.coord;
                    s.connec = NegPoissMesh.connec;
                    obj.mesh = Mesh.create(s);
                otherwise
                    obj.mesh = HexaMesh(2,1,1,20,5,5);
            end
        end

        function createMaterial(obj)
            mu     = ConstantFunction.create(1,obj.mesh);
            lambda = ConstantFunction.create(1,obj.mesh);
            ndim   = obj.mesh.ndim;
            kappa = IsotropicElasticMaterial.computeKappaFromShearAndLambda(mu,lambda,ndim);

            s.type  = 'ISOTROPIC';
            s.ndim  = ndim;
            s.bulk  = kappa;
            s.shear = mu;
            obj.material = Material.create(s);
            obj.matProp.mu = mu;
            obj.matProp.lambda = lambda;
        end

        function createDisplacementFun(obj)
            obj.uFun = LagrangianFunction.create(obj.mesh, obj.mesh.ndim, 'P1');
            obj.uFun.setFValues(obj.uFun.fValues + 0);
        end

        function createFunctionals(obj)
            obj.createLinearElasticityFunctional();
            obj.createNeohookeanFunctional();
        end

        function createLinearElasticityFunctional(obj)
            s.material = obj.material;
            s.mesh     = obj.mesh;
            elas = LinearElasticityFunctional(s);
            obj.linearElasticityFun = elas;
        end

        function createNeohookeanFunctional(obj)
            s.material = obj.matProp;
            s.mesh     = obj.mesh;
            neo = NeohookeanFunctional(s);
            obj.neohookeanFun = neo;
        end

        function setMonitoring(obj,cParams)
            s.shallDisplay = true;
            s.shallPrint   = true;
            obj.monitor = HyperelasticityMonitoring(s);
        end

        function createBoundaryConditions(obj)
            nSteps = 25;
            maxVal = 1;
            cParams.bcValues = linspace(0,maxVal,nSteps);
            cParams.type = 'ForceTractionX';
            bc = HyperelasticityBoundaryCreator(obj.mesh,cParams);
            obj.boundaryConditions = bc;
        end

        function [u,bc] = updateBoundaryConditions(obj,u)
            bc = obj.boundaryConditions.nextStep();
            u.setFValues(obj.updateInitialDisplacement(u,bc));
        end

        function u = updateInitialDisplacement(~,uOld,bc)
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







        function printFile(obj,iStep, u, reac)
            if obj.print
                fun = {u, reac};
                funNames = {'Displacement', 'Reactions'};
                a.mesh     = obj.mesh;
                a.filename = ['SIM_',obj.fileNameOut,'_',int2str(iStep)];
                a.fun      = fun;
                a.funNames = funNames;
                a.type     = 'Paraview';
                pst = FunctionPrinter.create(a);
                pst.print();
            end
        end

        function reacFun = computeReactions(obj,KR,u)
            uVal = obj.reshapeToVector(u.fValues);
            reactions = zeros(size(uVal));
            reactions(obj.constrainedDofs,1) = -KR*uVal;
            s.fValues = obj.reshapeToMatrix(reactions);
            s.order = 'P1';
            s.mesh  = obj.mesh;
            reacFun = LagrangianFunction(s);
        end
        % 
        % function plotStep(obj,intEnergy,nIterPerStep,iStep, normDisp, normReac)
        % 
        %     f = figure(1);
        %     clf(f)
        %     subplot(1,3,1)
        %     plot(intEnergy)
        %     title('Energies')
        %     xlabel('iteration')
        %     ylabel('Energy')
        %     subplot(1,3,2)
        %     plot(normDisp, normReac)
        %     title('Displacement-reaction')
        %     xlabel('displacement')
        %     ylabel('reactions')
        %     hold on
        %     subplot(1,3,3)
        %     bar(1:iStep, nIterPerStep)
        %     title('Number of iterations to converge $\Delta$ F')
        %     xlabel('step')
        %     ylabel('num. iterations per step')
        %     hold on
        %     drawnow
        % end

        function hasNot = hasNotConverged(obj,energy)
            TOL = 10e-12;
            energyOld = energy(end-1);
            energy    = energy(end);
            hasNot = abs(energy-energyOld)/energy > TOL;
        end



        function u = computeInitialElasticGuess(obj, bc)
            s.mesh = obj.mesh;
            s.scale = 'MACRO';
            s.material = obj.material;
            s.dim = '2D';
            s.boundaryConditions = bc;
            s.solverType = 'REDUCED';
            s.solverMode = 'DISP';
            s.solverCase = 'DIRECT';
            fem = ElasticProblem(s);
            fem.solve();
            %             u = obj.reshapeToVector(fem.uFun.fValues);
            u = fem.uFun;
        end





        function [Hff,Hr] = computeHessian(obj,u)
            H = obj.linearElasticityFun.computeHessian(u);
            %H = obj.neohookeanFun.computeHessian(u);
            f = obj.freeDofs;
            r = obj.constrainedDofs;
            Hff = H(f,f);
            Hr  = H(r,:);
        end

        function Rf = computeResidual(obj,u)
            Fext          = obj.computeExternalForces();
            [Fint,FintEl] = obj.computeInternalForces(u);
            R = Fint - Fext;% - react;
            f = obj.freeDofs;
            Rf = R(f);
        end

        function nIter = computeNumberOfIterations(obj,iter,iterOld,iStep)
            if iStep == 1
                nIter = iter-1;
            else
                nIter = iter-1-sum(iterOld(1:(iStep-1)));
            end
        end

        function computeFreeDofs(obj)
            bc = obj.createBoundaryConditions(1);
            nDofs = obj.uFun.nDofs;
            obj.freeDofs = setdiff(1:nDofs,bc.dirichlet_dofs);
            obj.constrainedDofs = bc.dirichlet_dofs;
        end



        function u = applyDirichletToUFun(obj,u,bc)
            u_k = reshape(u.fValues',[u.nDofs,1]);
            u_k(bc.dirichlet_dofs) = bc.dirichlet_vals;
            u.setFValues(reshape(u_k,[obj.mesh.ndim,obj.mesh.nnodes])');
        end

        function [Fext] = computeExternalForces(obj,perc)
            %             Fext = perc*obj.FextInitial;
            %             Fext = obj.reshapeToVector(Fext);
            Fext = zeros(obj.uFun.nDofs,1);
        end

        function [intFor,intForel] = computeInternalForces(obj,u)
            %intFor = obj.neohookeanFun.computeGradient(u);
            intFor = obj.linearElasticityFun.computeGradient(u);
        end

        function rshp = reshapeToVector(obj, A)
            rshp = reshape(A',[obj.uFun.nDofs,1]);
        end

        function rshp = reshapeToMatrix(obj, A)
            rshp = reshape(A,[obj.mesh.ndim,obj.mesh.nnodes])';
        end

        function plotVector(obj, vect)
            s.fValues = obj.reshapeToMatrix(vect);
            s.mesh = obj.mesh;
            s.order = 'P1';
            a = LagrangianFunction(s);
            a.plot;
        end


        % Other useful stuff later on
        function lambdas = computeStretches(obj)
            GradU2 = ActualGrad(obj.uFun);
            Id = Identity(obj.uFun);
            F = GradU2 + Id;
            C = F'*F;
            lambdas = (Eigen(C).^0.5);
        end

        function sigma = computeCauchyStress(obj)
            mu = obj.material.mu;
            lambda = obj.material.lambda;
            GradU2 = ActualGrad(obj.uFun);
            Id = Identity(obj.uFun);
            F = GradU2 + Id;
            b = F*F';
            jac = Det(F);
            sigma = mu.*(b-Id)./jac + lambda.*(log(jac)).*Id./jac;
            sigma.ndimf = [2 2];
        end


    end

end
