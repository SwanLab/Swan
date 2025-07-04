classdef HyperelasticProblem < handle

    properties (Access = public)
        uFun
        rFun
    end

    properties (Access = private)
        mesh
        neohookeanFun, linearElasticityFun
        material, materialElastic
        bc_case
        fileName
        meshGen
        printing
        bcApplier
        freeDofs, constrainedDofs
    end

    methods (Access = public)

        function obj = HyperelasticProblem(cParams)
            close all;
            obj.init(cParams);
            obj.computeFreeDofs();

            u  = obj.uFun;
            f = animatedline;

            nsteps = cParams.nsteps;
            iter = 1;
            nIterPerStep = [];

            reactions = zeros(u.nDofs,1);
            for iStep = 1:nsteps
                loadPercent = iStep/nsteps;
                bc = obj.createBoundaryConditions(loadPercent);
                % u  = obj.computeInitialElasticGuess(bc);
                u  = obj.applyDirichletToUFun(u,bc);

                intEnergy(iter) = obj.neohookeanFun.compute(u);
                % intEnergy(iter) = obj.linearElasticityFun.compute(u);
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

        function printFile(obj,iStep, u, reac)
            if obj.printing
                fun = {u, reac};
                funNames = {'Displacement', 'Reactions'};
                a.mesh     = obj.mesh;
                a.filename = ['SIM_',obj.fileName,'_',int2str(iStep)];
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

        function plotStep(obj,intEnergy,nIterPerStep,iStep, normDisp, normReac)

            f = figure(1);
            clf(f)
            subplot(1,3,1)
            plot(intEnergy)
            title('Energies')
            xlabel('iteration')
            ylabel('Energy')
            subplot(1,3,2)
            plot(normDisp, normReac)
            title('Displacement-reaction')
            xlabel('displacement')
            ylabel('reactions')
            hold on
            subplot(1,3,3)
            bar(1:iStep, nIterPerStep)
            title('Number of iterations to converge $\Delta$ F')
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

        function init(obj,cParams)
            obj.printing = cParams.printing;
            obj.bc_case  = cParams.bcCase;
            obj.fileName = cParams.fileName;
            obj.meshGen  = cParams.meshGen;
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

        function [Hff,Hr] = computeHessian(obj,u)
            H = obj.neohookeanFun.computeHessian(u);
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

        function createMesh(obj)
            switch obj.meshGen
                case {'Hole', 'HoleDirich'}
                    IM = Mesh.createFromGiD('hole_mesh_quad.m');
                    obj.mesh = IM;
                case {'Bending', 'Traction'}
                    obj.mesh = UnitQuadMesh(20,20);
                case {'Metamaterial'}
                    load('NegPoissMesh.mat')
                    s.coord = NegPoissMesh.coord;
                    s.connec = NegPoissMesh.connec;
                    obj.mesh = Mesh.create(s);
                case 'EIFEMMesh'
                    load(obj.fileName)
                    s.coord  = EIFEoper.MESH.COOR;          
                    s.connec = EIFEoper.MESH.CN;
                    mS       = Mesh.create(s);
                    obj.mesh = mS;
                otherwise
                    obj.mesh = HexaMesh(2,1,1,20,5,5);
                    % obj.mesh = UnitHexaMesh(15,15,15);
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
            E1  = Isotropic2dElasticMaterial.computeYoungFromShearAndBulk(G,K,N)  ;
            nu1 = Isotropic2dElasticMaterial.computePoissonFromFromShearAndBulk(G,K,N);
            E2        = G*(3*L+2*G)/(L+G);
            nu2       = L / (2*(L+G));
            E         = ConstantFunction.create(E1,obj.mesh);
            nu        = ConstantFunction.create(nu1,obj.mesh);
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
            obj.uFun.setFValues(obj.uFun.fValues + 0);
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


    end

end
