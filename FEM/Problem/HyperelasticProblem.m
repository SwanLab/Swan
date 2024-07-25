classdef HyperelasticProblem < handle
    
    properties (Access = public)
        uFun
    end

    properties (Access = private)
        mesh
        neohookeanFun, linearElasticityFun
        material, materialElastic
        bc_case = 'Traction'
    end

    methods (Access = public)

        function obj = HyperelasticProblem()
            close all;
            obj.init();

            % Apply boundary conditions
            bc = obj.createBoundaryConditions(1);
            nDofs = obj.uFun.nDofs;
            freeDofs = setdiff(1:nDofs,bc.dirichlet_dofs);
            u_k = reshape(obj.uFun.fValues',[obj.uFun.nDofs,1]);
            u_k(bc.dirichlet_dofs) = bc.dirichlet_vals;

            % Init Newton-Raphson
            f = animatedline;
%             fig2 = animatedline;
            obj.applyDirichletToUFun(bc);

            nsteps = 100;
            energies = [];
            alpha = 1e-3;
            react = zeros(obj.uFun.nDofs,1);
            for iStep = 1:nsteps
                disp('------')
                disp('NEW LOAD STEP')

                loadPercent = iStep/nsteps;
                bc = obj.createBoundaryConditions(loadPercent);
                obj.applyDirichletToUFun(bc);
                u_k = reshape(obj.uFun.fValues',[obj.uFun.nDofs,1]);
                Fext = obj.computeExternalForces();
                Fint = obj.computeInternalForces();
                % hess = neo.computeHessian(obj.uFun);
                R = Fint - Fext - react;

                u_k(freeDofs) = u_k(freeDofs) - alpha*R(freeDofs);
                nrg0 = obj.neohookeanFun.compute(obj.uFun);
                old_nrg = nrg0;

                residual = norm(R(freeDofs));
                resi = 0;
                i = 0;
                prenorm = 100;
                hasNotConverged = 1;
                while hasNotConverged %residual > 10e-8
                    % Update U
                    % [deltaUk,~] = obj.solveProblem(hess,-R,u_k);
                    
                    u_next = u_k + alpha*R;
%                     react = react + deltaReactk;
                    obj.uFun.fValues = obj.reshapeToMatrix(u_next);
                    obj.applyDirichletToUFun(bc);
                    u_k = obj.reshapeToVector(obj.uFun.fValues);
%                     u_k = u_next;
                    if ~isreal(obj.uFun.fValues)
                        warning('Complex numbers')
                    end

                    % Calculate residual
                    Fint = obj.computeInternalForces();
%                     hess = obj.neohookeanFun.computeHessian(obj.uFun);
                    R = Fint - Fext ;
                    nrg = obj.neohookeanFun.compute(obj.uFun);
                    hasNotConverged = abs(nrg-old_nrg)/nrg0 > 10e-12;
                    old_nrg = nrg;

                    % Plot
                    i = i+1;
                    abs(prenorm-residual)
                    resi(i) = residual;
                    energies(i) = nrg;
                    prenorm = residual;
                    
                end
                obj.uFun.print(['SIM_Bending_',int2str(iStep)])

                num_is(iStep) = i;
                f = figure(1);
                clf(f)
                subplot(1,2,1)
                plot(energies)
                title('Energies')
                xlabel('displacement (m)')
                ylabel('force (N)')
                subplot(1,2,2)
                bar(1:iStep, num_is)
                title('Number of iterations to converge \DeltaF')
                xlabel('step')
                ylabel('num. iterations')
                hold on
                drawnow
            end

        end

    end

    methods (Access = private)

        function init(obj)
            obj.createMesh();
            obj.createMaterial();
            obj.createDisplacementFun();
            obj.createFunctionals();
        end

        function fem = solveElasticProblem(obj, perc)
            s.mesh = obj.mesh;
            s.scale = 'MACRO';
            s.material = obj.materialElastic;
            s.dim = '2D';
            s.boundaryConditions = obj.createBoundaryConditions(perc);
            s.solverType = 'REDUCED';
            s.solverMode = 'DISP';
            fem = ElasticProblem(s);
            fem.solve();
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
        
        function createMesh(obj)
            switch obj.bc_case
                case {'Hole', 'HoleDirich'}
                    IM = Mesh.createFromGiD('hole_mesh_quad.m');
                    obj.mesh = IM;
                case {'Bending', 'Traction'}
                    obj.mesh = UnitQuadMesh(14,14);
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
            E1        = G*(3*L+2*G)/(L+G);
            nu1       = L / (2*(L+G));
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

        function applyDirichletToUFun(obj, bc)
            u_k = reshape(obj.uFun.fValues',[obj.uFun.nDofs,1]);
            u_k(bc.dirichlet_dofs) = bc.dirichlet_vals;
            obj.uFun.fValues = reshape(u_k,[obj.mesh.ndim,obj.mesh.nnodes])';
        end

        function [Fext] = computeExternalForces(obj,perc)
%             Fext = perc*obj.FextInitial;
%             Fext = obj.reshapeToVector(Fext);
            Fext = zeros(obj.uFun.nDofs,1);
        end

        function intfor = computeInternalForces(obj)
            intfor = obj.neohookeanFun.computeGradient(obj.uFun);
            intforel = obj.linearElasticityFun.computeGradient(obj.uFun);
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
