classdef HyperelasticProblem < handle
    
    properties (Access = public)
        uFun
        strainFun
        stressFun
    end

    properties (Access = private)
        quadrature
        boundaryConditions, BCApplier
        dirichletFun

        neohookeanFun

        stiffness
        Fext
        solver, solverType, solverMode, solverCase
        scale
        
        strain, stress

        FextInitial
    end

    properties (Access = protected)
        mesh 
        material
        displacementFun
    end

    methods (Access = public)

        function obj = HyperelasticProblem()
            close all;
            obj.init();
            obj.createDisplacementFun();
            obj.createBoundaryConditions();

            % Create Neohookean Functional
            s.material = obj.material;
            s.mesh     = obj.mesh;
            neo = NeohookeanFunctional(s);
            obj.neohookeanFun = neo;

            % Apply boundary conditions
            bc = obj.boundaryConditions;
            u_k = reshape(obj.uFun.fValues',[obj.uFun.nDofs,1]);
            u_k(bc.dirichlet_dofs) = bc.dirichlet_vals;

            % Init Newton-Raphson
            f = animatedline;
%             fig2 = animatedline;
            obj.applyDirichletToUFun();

            nsteps = 100;
            displ_grafic = [];
            fext_grafic = [];
            nrg_grafic = [];
            react = zeros(obj.uFun.nDofs,1);
            hess = neo.computeHessian(obj.uFun);
            for iStep = 1:nsteps
                disp('------')
                disp('NEW LOAD STEP')

                loadPercent = iStep/nsteps;
                obj.createBoundaryConditions();
                Fext = obj.computeExternalForces(loadPercent);
                Fint = obj.computeInternalForces();
                hess = neo.computeHessian(obj.uFun);
                R = Fint - Fext - react;

                residual = norm(R);
                resi = 0;
                i = 0;
                while residual > 10e-8
                    nrg = neo.compute(obj.uFun)

                    % Update U
                    [deltaUk,react_k] = obj.solveProblem(hess,-R);
                    
                    u_next = u_k + deltaUk;
                    react = react + react_k;
                    obj.uFun.fValues = reshape(u_next,[obj.mesh.ndim,obj.mesh.nnodes])';
                    u_k = u_next;
                    obj.uFun.fValues

                    % Calculate residual
                    Fint = obj.computeInternalForces();
                    hess = neo.computeHessian(obj.uFun);
                    R = Fint - Fext - react;

                    residual = norm(R)/norm(Fext)
                    sigma = obj.computeCauchyStress();
%                     lambdas = obj.computeStretches();

                    % Plot
                    i = i+1;
                    resi(i) = residual;
%                     addpoints(f,i,log10(residual));
%                     drawnow
                    
                end
                obj.uFun.print(['AAShapeFine_paraview',int2str(iStep)])

                xMax    = max(obj.mesh.coord(:,1));

                isRight  = @(coor)  abs(coor(:,1))==xMax;
                isBottom = @(coor)  abs(coor(:,2))==0;
                bot_right = @(coor) isBottom(coor) & isRight(coor);
                dof = obj.uFun.getDofsFromCondition(bot_right);
                dof = dof(1);
                nDimf = obj.uFun.ndimf;
                nNode = obj.mesh.nnodes;
                dofToDim = repmat(1:nDimf,[1,nNode]);
                dofToNode = repmat(1:nNode, nDimf, 1);
                dofToNode = dofToNode(:);

                nod = dofToNode(dof);
                dof_dir = dofToDim(dof); 

                nrg_grafic  = [nrg_grafic, nrg];
                displ_grafic = [displ_grafic, obj.uFun.fValues(nod, dof_dir)];
                fext_grafic  = [fext_grafic, Fext(dof)];
                posgp = Quadrature.create(obj.mesh,1).posgp;

                num_is(iStep) = i;
                f = figure(1);
                clf(f)
                subplot(1,3,1)
                plot(displ_grafic, fext_grafic,'-x')
                subplot(1,3,2)
                plot(1:iStep, nrg_grafic,'-x')
                subplot(1,3,3)
                bar(1:iStep, num_is)
                hold on
                drawnow
%                 addpoints(fig2, obj.uFun.fValues(3,1), Fext(5))
%                 drawnow
%                 num_is(iStep) = i;

            end

        end

        function createBCApplier(obj)
            s.mesh = obj.mesh;
            s.boundaryConditions = obj.boundaryConditions;
            bc = BCApplier(s);
            obj.BCApplier = bc;
        end

        function [u,reactions] = solveProblem(obj, lhs, rhs)
            obj.createBCApplier();
            a.type = 'DIRECT';
            solv = Solver.create(a);
            s.solverType = 'MONOLITHIC';
            s.solverMode = 'DISP';
            s.stiffness  = lhs;
            s.forces     = rhs;
            s.solver     = solv;
            s.boundaryConditions = obj.boundaryConditions;
            s.BCApplier          = obj.BCApplier;
            pb = ProblemSolver(s);
            [u,L] = pb.solve();

            reactions = zeros(obj.uFun.nDofs, 1);
            reactions(obj.boundaryConditions.dirichlet_dofs) = L;
            reac_rshp = reshape(reactions,[obj.mesh.ndim,obj.mesh.nnodes])';
        end

    end

    methods (Access = private)

        function init(obj)
%             obj.mesh = HexaMesh(2,1,1,20,5,5);
%             obj.mesh = UnitHexaMesh(5,5,5);
            obj.mesh = UnitQuadMesh(50,50);

%             obj.material.mu = 3/8;
%             obj.material.lambda = 3/4;

            obj.material.mu = 1;
            obj.material.lambda = 1*1.8;

            k = obj.material.lambda + 2/obj.mesh.ndim * obj.material.mu;
            m = obj.material.mu;

            N = obj.mesh.ndim;
            E = ((N*N*k).*(2*m))./(2*m + N*(N-1)*k);
            nu = ((N*k)-(2*m))./(2*m + N*(N-1)*k);

%             E = 10.0;
%             nu = 0.3;
% 
            mu = E/(2*(1 + nu));
            k = E./(N*(1-(N-1)*nu));
            lambda = k - 2/N*mu;

            obj.material.lambda = lambda;
            obj.material.mu = mu;
        end
        
        function createBoundaryConditions(obj)
%             obj.createBC2D_oneelem();
%             obj.createBC2D_knownexample();
%             obj.createBCflexio();
            obj.createBC2D_nelem();
%             obj.createBC3D_oneelem();
%             obj.createBC3D_nelem();
        end

        function createDisplacementFun(obj)
            obj.uFun = LagrangianFunction.create(obj.mesh, obj.mesh.ndim, 'P1');
            obj.uFun.fValues = obj.uFun.fValues + 0;
        end

        function applyDirichletToUFun(obj)
            bc = obj.boundaryConditions;
            u_k = reshape(obj.uFun.fValues',[obj.uFun.nDofs,1]);
            u_k(bc.dirichlet_dofs) = bc.dirichlet_vals;
            obj.uFun.fValues = reshape(u_k,[obj.mesh.ndim,obj.mesh.nnodes])';
        end

        function [Fext] = computeExternalForces(obj,perc)
%             pl = obj.boundaryConditions.pointloadFun;
%             pl.fValues = pl.fValues*perc;
%             s.mesh = obj.mesh;
%             s.type = 'ShapeFunction';
%             s.quadType = 2;
%             rhsI       = RHSintegrator.create(s);
%             test = LagrangianFunction.create(obj.mesh,obj.uFun.ndimf,'P1');
%             Fext2 = rhsI.compute(pl,test);
% 
%             Fext = pl.fValues;
%             Fext = reshape(Fext',[obj.uFun.nDofs,1]);
% %             sum(Fext)
% 
%             s.mesh = obj.mesh;
%             s.quadType = 2;
%             scalar = IntegratorScalarProduct(s);
%             intF = scalar.compute(pl,pl)
% %             Fext = reshape(Fext2,[obj.mesh.ndim,obj.mesh.nnodes])';
%             sumFext2 = sum(Fext2)
%             Fext = Fext2;
            Fext = perc*obj.FextInitial;
            Fext = obj.reshapeToVector(Fext);
        end

        function intfor = computeInternalForces(obj)
            intfor = obj.neohookeanFun.computeInternalForces(obj.uFun,obj.boundaryConditions);
        end

        function lambdas = computeStretches(obj)
            GradU2 = ActualGrad(obj.uFun);
            Id = Identity(obj.uFun);
            F = GradU2 + Id;
            C = F'*F;
            lambdas = (Eigen(C).^0.5);
        end

        function rshp = reshapeToVector(obj, A)
            rshp = reshape(A',[obj.uFun.nDofs,1]);
        end

        function rshp = reshapeToMatrix(obj, A)
            rshp = reshape(A,[obj.mesh.ndim,obj.mesh.nnodes])';
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

%         function bc = createBC2D_oneelem(obj)
%             xMax    = max(obj.mesh.coord(:,1));
%             yMax    = max(obj.mesh.coord(:,2));
%             isLeft   = @(coor)  abs(coor(:,1))==0;
%             isRight  = @(coor)  abs(coor(:,1))==xMax;
%             isTop    = @(coor)  abs(coor(:,2))==yMax;
%             isBottom = @(coor)  abs(coor(:,2))==0;
%             isMiddle = @(coor)  abs(coor(:,2))==yMax/2;
% 
%             % 2D ONE ELEMENT
%             sDir1.domain    = @(coor) isLeft(coor) & isTop(coor);
%             sDir1.direction = [1];
%             sDir1.value     = 0;
%             dir1 =  DirichletCondition(obj.mesh, sDir1);
% 
%             sDir2.domain    = @(coor) isLeft(coor) & isBottom(coor);
%             sDir2.direction = [1,2];
%             sDir2.value     = 0;
%             dir2 =  DirichletCondition(obj.mesh, sDir2);
%             s.dirichletFun = [dir1, dir2];
% 
%             sPL.domain    = @(coor) isRight(coor);
%             sPL.direction = 1;
%             sPL.value     = 0.1;
%             s.pointloadFun = PointLoad(obj.mesh, sPL);
%             
%             s.periodicFun  = [];
%             s.mesh         = obj.mesh;
% 
%             bc = BoundaryConditions(s);
%             obj.boundaryConditions = bc;
%         end

%         function bc = createBC2D_knownexample(obj)
%             xMax    = max(obj.mesh.coord(:,1));
%             yMax    = max(obj.mesh.coord(:,2));
%             isLeft   = @(coor)  abs(coor(:,1))==0;
%             isRight  = @(coor)  abs(coor(:,1))==xMax;
%             isTop    = @(coor)  abs(coor(:,2))==yMax;
%             isBottom = @(coor)  abs(coor(:,2))==0;
%             isMiddle = @(coor)  abs(coor(:,2))==yMax/2;
% 
%             constant = 1.08*obj.material.lambda*0.5;
% 
%             % Node 1
%             sDir1xy.domain    = @(coor) isLeft(coor) & isBottom(coor);
%             sDir1xy.direction = [1,2];
%             sDir1xy.value     = 0;
%             dir1 =  DirichletCondition(obj.mesh, sDir1xy);
% 
%             % Right, x
%             sDir2x.domain    = @(coor) isRight(coor);
%             sDir2x.direction = [1];
%             sDir2x.value     = 1.2;
%             dir2x =  DirichletCondition(obj.mesh, sDir2x);
% 
%             % Bottom right, y
%             sDir2a.domain    = @(coor) isBottom(coor);
%             sDir2a.direction = [2];
%             sDir2a.value     = 0;
%             dir2aa =  DirichletCondition(obj.mesh, sDir2a);
% 
%             % Top, y
%             sDir2y.domain    = @(coor) isTop(coor);
%             sDir2y.direction = [2];
%             sDir2y.value     = 0.9;
%             dir2y =  DirichletCondition(obj.mesh, sDir2y);
% 
%             % Top left, x
%             sDir2z.domain    = @(coor) isTop(coor) & isLeft(coor);
%             sDir2z.direction = [1];
%             sDir2z.value     = 0;
%             dir2z =  DirichletCondition(obj.mesh, sDir2z);
% 
%             s.dirichletFun = [dir1, dir2x, dir2aa, dir2y, dir2z];
% 
%             sPL.domain    = @(coor) isRight(coor);
%             sPL.direction = 1;
%             sPL.value     = 0;
%             s.pointloadFun = PointLoad(obj.mesh, sPL);
%             
%             s.periodicFun  = [];
%             s.mesh         = obj.mesh;
% 
%             bc = BoundaryConditions(s);
%             obj.boundaryConditions = bc;
%         end

%         function bc = createBCflexio(obj)
%             xMax    = max(obj.mesh.coord(:,1));
%             yMax    = max(obj.mesh.coord(:,2));
%             isLeft   = @(coor)  abs(coor(:,1))==0;
%             isRight  = @(coor)  abs(coor(:,1))==xMax;
%             isHalf   = @(coor)  abs(coor(:,1))==xMax/2;
%             isTop    = @(coor)  abs(coor(:,2))==yMax;
%             isBottom = @(coor)  abs(coor(:,2))==0;
%             isMiddle = @(coor)  abs(coor(:,2))==yMax/2;
% 
%             % 2D BENDING
%             sDir1.domain    = @(coor) isLeft(coor) & ~isBottom(coor);
%             sDir1.direction = [1];
%             sDir1.value     = 0;
%             dir1 =  DirichletCondition(obj.mesh, sDir1);
% 
%             sDir2.domain    = @(coor) isLeft(coor) & isBottom(coor);
%             sDir2.direction = [1,2];
%             sDir2.value     = 0;
%             dir2 =  DirichletCondition(obj.mesh, sDir2);
% 
%             sDir3.domain    = @(coor) isRight(coor) & ~isBottom(coor);
%             sDir3.direction = [1];
%             sDir3.value     = 0;
%             dir3 =  DirichletCondition(obj.mesh, sDir3);
% 
%             sDir4.domain    = @(coor) isRight(coor) & isBottom(coor);
%             sDir4.direction = [1,2];
%             sDir4.value     = 0;
%             dir4 =  DirichletCondition(obj.mesh, sDir4);
%             s.dirichletFun = [dir1, dir2, dir3, dir4];
% 
%             sPL.domain    = @(coor) isTop(coor) & isHalf(coor);
%             sPL.direction = 2;
%             sPL.value     = -0.1;
%             s.pointloadFun = PointLoad(obj.mesh, sPL);
%             
%             s.periodicFun  = [];
%             s.mesh         = obj.mesh;
% 
%             bc = BoundaryConditions(s);
%             obj.boundaryConditions = bc;
%         end

        function bc = createBC2D_nelem(obj)
            xMax    = max(obj.mesh.coord(:,1));
            yMax    = max(obj.mesh.coord(:,2));
            isLeft   = @(coor)  abs(coor(:,1))==0;
            isRight  = @(coor)  abs(coor(:,1))==xMax;
            isTop    = @(coor)  abs(coor(:,2))==yMax;
            isBottom = @(coor)  abs(coor(:,2))==0;
            isMiddle = @(coor)  abs(abs(coor(:,2))-yMax/2) <= 10e-2;
            
            % 2D N ELEMENTS
            sDir1.domain    = @(coor) isLeft(coor) & ~isMiddle(coor);
            sDir1.direction = [1];
            sDir1.value     = 0;
            dir1 =  DirichletCondition(obj.mesh, sDir1);

            sDir2.domain    = @(coor) isLeft(coor) & isMiddle(coor);
            sDir2.direction = [1,2];
            sDir2.value     = 0;
            dir2 =  DirichletCondition(obj.mesh, sDir2);
            s.dirichletFun = [dir1, dir2];

            sPL.domain    = @(coor) isRight(coor);
            sPL.direction = 1;
            sPL.value     = 1;
            s.pointloadFun = [];DistributedLoad(obj.mesh, sPL);

            [bM,l2g] = obj.mesh.getBoundarySubmesh(sPL.domain);

            sAF.fHandle = @(x) [sPL.value*ones(size(x(1,:,:)));0*x(2,:,:)];
            sAF.ndimf   = 2;
            sAF.mesh    = bM;
            xFun = AnalyticalFunction(sAF);
            xFunP1  =xFun.project('P1');

            s.mesh = bM;
            s.type = 'ShapeFunction';
            s.quadType = 2;
            rhsI       = RHSintegrator.create(s);
            test = LagrangianFunction.create(bM,xFun.ndimf,'P1');
            Fext2 = rhsI.compute(xFunP1,test);   
            Fext3 = reshape(Fext2,[bM.ndim,bM.nnodes])';

            Fext = zeros(obj.mesh.nnodes,2);
            Fext(l2g,:) = Fext3;

            obj.FextInitial = Fext; 
            
            s.periodicFun  = [];
            s.mesh         = obj.mesh;

            bc = BoundaryConditions(s);
            obj.boundaryConditions = bc;
        end

%         function bc = createBC3D_oneelem(obj)
%             xMax    = max(obj.mesh.coord(:,1));
%             zMax    = max(obj.mesh.coord(:,3));
%             isLeft   = @(coor)  abs(coor(:,1))==0;
%             isRight  = @(coor)  abs(coor(:,1))==xMax;
%             isTop    = @(coor)  abs(coor(:,3))==zMax;
%             isBottom = @(coor)  abs(coor(:,3))==0;
%             isMiddle = @(coor)  abs(coor(:,3))==zMax/2;
% 
%             % 2D ONE ELEMENT
%             sDir1.domain    = @(coor) isLeft(coor) & isTop(coor);
%             sDir1.direction = [1];
%             sDir1.value     = 0;
%             dir1 =  DirichletCondition(obj.mesh, sDir1);
% 
%             sDir2.domain    = @(coor) isLeft(coor) & isBottom(coor);
%             sDir2.direction = [1,2,3];
%             sDir2.value     = 0;
%             dir2 =  DirichletCondition(obj.mesh, sDir2);
%             s.dirichletFun = [dir1, dir2];
% 
%             sPL.domain    = @(coor) isRight(coor);
%             sPL.direction = 1;
%             sPL.value     = 0.1;
%             s.pointloadFun = PointLoad(obj.mesh, sPL);
%             
%             s.periodicFun  = [];
%             s.mesh         = obj.mesh;
% 
%             bc = BoundaryConditions(s);
%             obj.boundaryConditions = bc;
%         end
% 
%         function bc = createBC3D_nelem(obj)
%             xMax    = max(obj.mesh.coord(:,1));
%             zMax    = max(obj.mesh.coord(:,3));
%             isLeft   = @(coor)  abs(coor(:,1))==0;
%             isRight  = @(coor)  abs(coor(:,1))==xMax;
%             isTop    = @(coor)  abs(coor(:,3))==zMax;
%             isBottom = @(coor)  abs(coor(:,3))==0;
%             isMiddle = @(coor)  abs(coor(:,3))==zMax/2;
%             
%             % 3D N ELEMENTS
% %             sDir1.domain    = @(coor) isLeft(coor) & ~isMiddle(coor);
% %             sDir1.direction = [1];
% %             sDir1.value     = 0;
% %             dir1 =  DirichletCondition(obj.mesh, sDir1);
% % 
% %             sDir2.domain    = @(coor) isLeft(coor) & isMiddle(coor);
% %             sDir2.direction = [1,2,3];
% %             sDir2.value     = 0;
% %             dir2 =  DirichletCondition(obj.mesh, sDir2);
% %             s.dirichletFun = [dir1, dir2];
% 
%             sDir2.domain    = @(coor) isLeft(coor);
%             sDir2.direction = [1,2,3];
%             sDir2.value     = 0;
%             s.dirichletFun =  DirichletCondition(obj.mesh, sDir2);
% 
%             sPL.domain    = @(coor) isRight(coor);
%             sPL.direction = 1;
%             sPL.value     = 10;
%             s.pointloadFun = PointLoad(obj.mesh, sPL);
%             
%             s.periodicFun  = [];
%             s.mesh         = obj.mesh;
% 
%             bc = BoundaryConditions(s);
%             obj.boundaryConditions = bc;
%         end

    end

end
