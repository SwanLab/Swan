classdef HyperelasticProblem < handle

    properties (Access = public)
        uFun
    end

    properties (Access = private)
        mesh
        neohookeanFun, linearElasticityFun
        material, materialElastic
        bc_case = 'HoleDirich'
        bcApplier
        freeDofs
    end

    methods (Access = public)



        function obj = HyperelasticProblem()
            close all;
            obj.init();
            obj.computeFreeDofs();

            u  = obj.uFun;
            f = animatedline;

            nsteps = 50;
            iter = 1;
            nIterPerStep = [];

            for iStep = 1:nsteps
                disp('------')
                disp('NEW LOAD STEP')

                loadPercent = iStep/nsteps;
                bc = obj.createBoundaryConditions(loadPercent);
                u  = obj.computeInitialElasticGuess(bc);
                u  = obj.applyDirichletToUFun(u,bc);

%                 intEnergy(iter) = obj.neohookeanFun.compute(obj.uFun);
                intEnergy(iter) = obj.linearElasticityFun.compute(obj.uFun);
                Res = obj.computeResidual(u);
                resi(iter) = norm(Res);



                hasNotConverged = 1;
                while hasNotConverged
                    % Update U

                    uVal = obj.reshapeToVector(u.fValues);
                    hess = obj.computeHessian(u);
                    Res  = obj.computeResidual(u);
                    incU = hess\(-Res);
                    uVal(obj.freeDofs) = uVal(obj.freeDofs) + incU;
                    u.fValues = obj.reshapeToMatrix(uVal);


                    iter = iter+1;
                    intEnergy(iter) = obj.neohookeanFun.compute(u);
                    Res        = obj.computeResidual(u);
                    resi(iter) = norm(Res);

                    hasNotConverged = obj.hasNotConverged(intEnergy);


                end
                obj.uFun.print(['SIM_',obj.bc_case,'_',int2str(iStep)])                
                nIterPerStep(iStep) = obj.computeNumberOfIterations(iter,nIterPerStep,iStep);
                obj.plotStep(intEnergy,nIterPerStep,iStep);
            end

        end

    end

    methods (Access = private)

        function plotStep(obj,intEnergy,nIterPerStep,iStep)

            f = figure(1);
            clf(f)
            subplot(1,2,1)
            plot(intEnergy)
            title('Energies')
            xlabel('iteration (m)')
            ylabel('Energy (N)')
            subplot(1,2,2)
            bar(1:iStep, nIterPerStep)
            title('Number of iterations to converge \DeltaF')
            xlabel('step')
            ylabel('num. iterations per step')
            hold on
            drawnow
        end

        function hasNot = hasNotConverged(obj,energy)
            TOL = 10e-12;
            energyOld = energy(end-1);
            energy    = energy(end);
            hasNot = abs(energy-energyOld)/energy > TOL;
        end

        function init(obj)
            obj.createMesh();
            obj.createMaterial();
            obj.createDisplacementFun();
            obj.createFunctionals();
        end

        function u = computeInitialElasticGuess(obj, bc)
            s.mesh = obj.mesh;
            s.scale = 'MACRO';
            s.material = obj.materialElastic;
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

        function createFunctionals(obj)
            obj.createLinearElasticityFunctional();
            obj.createNeohookeanFunctional();
        end

        function createLinearElasticityFunctional(obj)
            s.material = obj.materialElastic;
            s.mesh     = obj.mesh;
            elas = LinearElasticityFunctional(s);
            obj.linearElasticityFun = elas;
        end

        function createNeohookeanFunctional(obj)
            s.material = obj.material;
            s.mesh     = obj.mesh;
            neo = NeohookeanFunctional(s);
            obj.neohookeanFun = neo;
        end

        function Hf = computeHessian(obj,u)
            H = obj.neohookeanFun.computeHessian(u);
            f = obj.freeDofs;
            Hf = H(f,f);
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
                nIter = iter-1-sum(iterOld(1:iStep-1));
            end
        end

        function computeFreeDofs(obj)
            bc = obj.createBoundaryConditions(1);
            nDofs = obj.uFun.nDofs;
            obj.freeDofs = setdiff(1:nDofs,bc.dirichlet_dofs);
        end

        function createMesh(obj)
            switch obj.bc_case
                case {'Hole', 'HoleDirich'}
                    IM = Mesh.createFromGiD('hole_mesh_quad.m');
                    obj.mesh = IM;
                case {'Bending', 'Traction'}
                    obj.mesh = UnitQuadMesh(20,20);
                otherwise
                    obj.mesh = HexaMesh(2,1,1,20,5,5);
                    %                     obj.mesh = UnitHexaMesh(15,15,15);
            end
        end

        function createMaterial(obj)
            obj.material.mu = 1;
            obj.material.lambda = 1*1;
            obj.materialElastic = obj.createElasticMaterial();
        end

        function mat = createElasticMaterial(obj)

            G = obj.material.mu;
            L = obj.material.lambda;
            N = obj.mesh.ndim;
            K = 2/N*G + L;
            E1  = Isotropic2dElasticMaterial.computeYoungFromShearAndBulk(G,K,N);
            nu1 = Isotropic2dElasticMaterial.computePoissonFromFromShearAndBulk(G,K,N);
            E2        = G*(3*L+2*G)/(L+G);
            nu2       = L / (2*(L+G));
            E         = AnalyticalFunction.create(@(x) E1*ones(size(squeeze(x(1,:,:)))),1,obj.mesh);
            nu        = AnalyticalFunction.create(@(x) nu1*ones(size(squeeze(x(1,:,:)))),1,obj.mesh);
            s.pdim    = obj.mesh.ndim;
            s.nelem   = obj.mesh.nelem;
            s.mesh    = obj.mesh;
            s.young   = E;
            s.poisson = nu;
            s.type = 'ISOTROPIC';
            s.ndim = obj.mesh.ndim;
            mat = Material.create(s);
        end

        function bcs = createBoundaryConditions(obj,perc)
            bc = HyperelasticityTestBCs(obj.bc_case, obj.mesh, perc);
            bcs = bc.boundaryConditions;
        end

        function createDisplacementFun(obj)
            obj.uFun = LagrangianFunction.create(obj.mesh, obj.mesh.ndim, 'P1');
            obj.uFun.fValues = obj.uFun.fValues + 0;
        end

        function u = applyDirichletToUFun(obj,u,bc)
            u_k = reshape(u.fValues',[u.nDofs,1]);
            u_k(bc.dirichlet_dofs) = bc.dirichlet_vals;
            u.fValues = reshape(u_k,[obj.mesh.ndim,obj.mesh.nnodes])';
        end

        function [Fext] = computeExternalForces(obj,perc)
            %             Fext = perc*obj.FextInitial;
            %             Fext = obj.reshapeToVector(Fext);
            Fext = zeros(obj.uFun.nDofs,1);
        end

        function [intFor,intForel] = computeInternalForces(obj,u)
            intFor = obj.neohookeanFun.computeGradient(u);
            intForel = obj.linearElasticityFun.computeGradient(u);
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

        function [incU,incR] = solveProblem(obj, lhs, rhs,uOld,bc)
            obj.createBCApplier();
            a.type = 'DIRECT';
            solv = Solver.create(a);
            s.solverType = 'REDUCED';
            s.solverMode = 'DISP';
            s.stiffness  = lhs;
            s.forces     = rhs;
            s.uOld       = uOld;
            s.solver     = solv;
            s.boundaryConditions = bc;
            s.BCApplier          = obj.bcApplier;
            pb = ProblemSolver(s);
            [incU,incL] = pb.solve();

            incR = zeros(obj.uFun.nDofs, 1);
            incR(bc.dirichlet_dofs) =  incL;
            reac_rshp = reshape(incR,[obj.mesh.ndim,obj.mesh.nnodes])';
        end

        function createBCApplier(obj,bc)
            s.mesh = obj.mesh;
            s.boundaryConditions = bc;
            bc = BCApplier(s);
            obj.bcApplier = bc;
        end


    end

end
