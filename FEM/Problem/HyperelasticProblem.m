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
        %Fext
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
            obj.createBoundaryConditions(1);

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
            for iStep = 1:nsteps
                disp('------')
                disp('NEW LOAD STEP')

                loadPercent = iStep/nsteps;
                obj.createBoundaryConditions(loadPercent);
                obj.applyDirichletToUFun();
                u_k = reshape(obj.uFun.fValues',[obj.uFun.nDofs,1]);
                Fext = obj.computeExternalForces(loadPercent);
                Fint = obj.computeInternalForces();
                hess = neo.computeHessian(obj.uFun);
                R = Fint - Fext - react;

                residual = norm(R);
                resi = 0;
                i = 0;
                while residual > 10e-8
%                     nrg = neo.compute(obj.uFun)

                    % Update U
                    [deltaUk,react_k] = obj.solveProblem(hess,-R);
                    
                    u_next = u_k + deltaUk;
                    react = react + react_k;
                    obj.uFun.fValues = reshape(u_next,[obj.mesh.ndim,obj.mesh.nnodes])';
                    obj.applyDirichletToUFun();
                    u_k = obj.reshapeToVector(obj.uFun.fValues);
%                     u_k = u_next;
                    obj.uFun.fValues

                    % Calculate residual
                    Fint = obj.computeInternalForces();
                    hess = neo.computeHessian(obj.uFun);
                    R = Fint - Fext - react;

                    residual = norm(R)%/norm(Fext)
%                     sigma = obj.computeCauchyStress();
%                     lambdas = obj.computeStretches();

                    % Plot
                    i = i+1;
                    resi(i) = residual;
%                     addpoints(f,i,log10(residual));
%                     drawnow
                    
                end
                obj.uFun.print(['SIM_Bending_',int2str(iStep)])

                xMax    = max(obj.mesh.coord(:,1));
                yMax    = max(obj.mesh.coord(:,2));

                isRight  = @(coor)  abs(coor(:,1))==xMax;
                isBottom = @(coor)  abs(coor(:,2))==0;
                isLeft  = @(coor)  abs(coor(:,1))==0;
                isTop = @(coor)  abs(coor(:,2))==yMax;
                isMiddle  = @(coor)  abs(coor(:,1))==xMax/2;
%                 bot_right = @(coor) isBottom(coor) & isRight(coor);
                bot_right = @(coor) isMiddle(coor) & isTop(coor);
                dof = obj.uFun.getDofsFromCondition(bot_right);
                dof = dof(2);
                nDimf = obj.uFun.ndimf;
                nNode = obj.mesh.nnodes;
                dofToDim = repmat(1:nDimf,[1,nNode]);
                dofToNode = repmat(1:nNode, nDimf, 1);
                dofToNode = dofToNode(:);

                nod = dofToNode(dof);
                dof_dir = dofToDim(dof); 

                nrg_grafic  = [nrg_grafic, nrg];
                displ_grafic = [displ_grafic, obj.uFun.fValues(nod, dof_dir)];
                fext_grafic  = [fext_grafic, sum(Fext)];
                posgp = Quadrature.create(obj.mesh,1).posgp;

                num_is(iStep) = i;
                f = figure(1);
                clf(f)
                subplot(1,2,1)
                plot(displ_grafic, fext_grafic,'-x')
                title('Displacement of central top node')
                xlabel('displacement (m)')
                ylabel('force (N)')
%                 subplot(1,3,2)
%                 plot(1:iStep, nrg_grafic,'-x')
                subplot(1,2,2)
                bar(1:iStep, num_is)
                title('Number of iterations to converge \DeltaF')
                xlabel('step')
                ylabel('num. iterations')
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
            obj.createMesh();
            obj.createMaterial();
        end

        function createMesh(obj)
%             obj.mesh = HexaMesh(2,1,1,20,5,5);
%             obj.mesh = UnitHexaMesh(15,15,15);
%             obj.mesh = UnitQuadMesh(14,14);

            % Hole mesh
            mesh = UnitQuadMesh(20,20);
            gPar.type          = 'Circle';
            gPar.radius        = 0.25;
            gPar.xCoorCenter   = 0.5;
            gPar.yCoorCenter   = 0.5;
            g                  = GeometricalFunction(gPar);
            phiFun             = g.computeLevelSetFunction(mesh);
            lsCircle           = phiFun.fValues;
            lsCircleInclusion  = -lsCircle;
            sUm.backgroundMesh = mesh;
            sUm.boundaryMesh   = mesh.createBoundaryMesh;
            uMesh              = UnfittedMesh(sUm);
            uMesh.compute(lsCircleInclusion);
            
            IM = uMesh.createInnerMesh();
            obj.mesh = IM;
        end
        
        function createMaterial(obj)
            % Only 3D
            obj.material.mu = 1*1000;        % kPa
            obj.material.lambda = 1*10*1000; % kPa

            k = obj.material.lambda + 2/3 * obj.material.mu; % canviar
            G = obj.material.mu;
            L = obj.material.lambda;

            k = L + 2*G/3;
            E = G*(3*L+2*G)/(L+G);
            NU = L / (2*(L+G));

            % Using both 2D and 3D formulae
%             gg = E/(2*(1+NU));
%             kk = E/(3*(1-2*NU));
%             ll = E*NU/( (1+NU)*(1-2*NU));
% 
%             N = obj.mesh.ndim;
%             E = G*(3*obj.material.lambda+2*G)/(obj.material.lambda+G);
%             nu = obj.material.lambda/(2*(obj.material.lambda + k));
% 
% %             E = 10.0;
% %             nu = 0.3;
% % 
%             mu = E/(2*(1 + nu));
%             k = E./(N*(1-(N-1)*nu));
%             lambda = k - 2/N*mu;
% 
%             obj.material.lambda = lambda;
%             obj.material.mu = mu;
        end

        function createBoundaryConditions(obj,perc)
%             obj.createBC_2DTraction();
%             obj.createBC_2DBending();
%             obj.createBC_2DHole();
              obj.createBC_2DHoleDirich(perc);
%             obj.createBC_3DCube();
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

        function bc = createBC_2DTraction(obj)
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
            s.pointloadFun = [];%DistributedLoad(obj.mesh, sPL);

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

        function bc = createBC_2DBending(obj)
            xMax    = max(obj.mesh.coord(:,1));
            yMax    = max(obj.mesh.coord(:,2));
            isLeft   = @(coor)  abs(coor(:,1))==0;
            isRight  = @(coor)  abs(coor(:,1))==xMax;
            isHalf   = @(coor)  abs(abs(coor(:,1)) - xMax/2) <=10e-2;
            isTop    = @(coor)  abs(coor(:,2))==yMax;
            isBottom = @(coor)  abs(coor(:,2))==0;
            isMiddle = @(coor)  abs(abs(coor(:,2))-yMax/2) <= 10e-2;
            
            % 2D N ELEMENTS
            sDir2.domain    = @(coor) isLeft(coor);
            sDir2.direction = [1,2];
            sDir2.value     = 0;
            dir2 =  DirichletCondition(obj.mesh, sDir2);

            sDir4.domain    = @(coor) isRight(coor);
            sDir4.direction = [1,2];
            sDir4.value     = 0;
            dir4 =  DirichletCondition(obj.mesh, sDir4);
            s.dirichletFun = [dir2, dir4];
% 
            sPL.domain    = @(coor) isTop(coor) & isHalf(coor);
            sPL.direction = 2;
            sPL.value     = -1000;
            s.pointloadFun = [];%DistributedLoad(obj.mesh, sPL);

            [bM,l2g] = obj.mesh.getBoundarySubmesh(sPL.domain);

            sAF.fHandle = @(x) [0*x(1,:,:);sPL.value*ones(size(x(1,:,:)))];
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

        function bc = createBC_2DHoleDirich(obj, perc)
            xMax    = max(obj.mesh.coord(:,1));
            yMax    = max(obj.mesh.coord(:,2));
            isLeft   = @(coor)  abs(coor(:,1))==0;
            isRight  = @(coor)  abs(coor(:,1))==xMax;
            isHalf   = @(coor)  abs(abs(coor(:,1)) - xMax/2) <=10e-2;
            isTop    = @(coor)  abs(coor(:,2))==yMax;
            isBottom = @(coor)  abs(coor(:,2))==0;
            isMiddle = @(coor)  abs(abs(coor(:,2))-yMax/2) <= 10e-2;
            
            % 2D N ELEMENTS
            sDir.domain    = @(coor) isBottom(coor);
            sDir.direction = [1,2];
            sDir.value     = 0;
            dir =  DirichletCondition(obj.mesh, sDir);

            sDir2.domain    = @(coor) isTop(coor);
            sDir2.direction = [2];
            sDir2.value     = perc*1;
            dir2 =  DirichletCondition(obj.mesh, sDir2);

            sDir3.domain    = @(coor) isTop(coor);
            sDir3.direction = [1];
            sDir3.value     = 0;
            dir3 =  DirichletCondition(obj.mesh, sDir3);

            s.dirichletFun = [dir, dir2, dir3];
%             s.dirichletFun = [dir, dir2];
            s.pointloadFun = [];%DistributedLoad(obj.mesh, sPL);

            Fext = zeros(obj.mesh.nnodes,2);

            obj.FextInitial = Fext; 
            
            s.periodicFun  = [];
            s.mesh         = obj.mesh;

            bc = BoundaryConditions(s);
            obj.boundaryConditions = bc;
        end

        function bc = createBC_2DHole(obj)
            xMax    = max(obj.mesh.coord(:,1));
            yMax    = max(obj.mesh.coord(:,2));
            isLeft   = @(coor)  abs(coor(:,1))==0;
            isRight  = @(coor)  abs(coor(:,1))==xMax;
            isHalf   = @(coor)  abs(abs(coor(:,1)) - xMax/2) <=10e-2;
            isTop    = @(coor)  abs(coor(:,2))==yMax;
            isBottom = @(coor)  abs(coor(:,2))==0;
            isMiddle = @(coor)  abs(abs(coor(:,2))-yMax/2) <= 10e-2;
            
            % 2D N ELEMENTS
            sDir2.domain    = @(coor) isBottom(coor);
            sDir2.direction = [1,2];
            sDir2.value     = 0;
            dir2 =  DirichletCondition(obj.mesh, sDir2);
            s.dirichletFun = [dir2];
% 
            sPL.domain    = @(coor) isTop(coor);
            sPL.direction = 2;
            sPL.value     = 1;
            s.pointloadFun = [];%DistributedLoad(obj.mesh, sPL);

            [bM,l2g] = obj.mesh.getBoundarySubmesh(sPL.domain);

            sAF.fHandle = @(x) [0*x(1,:,:);sPL.value*ones(size(x(1,:,:)))];
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

        function bc = createBC_3DCube(obj)
            xMax    = max(obj.mesh.coord(:,1));
            yMax    = max(obj.mesh.coord(:,3));
            zMax    = max(obj.mesh.coord(:,3));
            isLeft   = @(coor)  abs(coor(:,1))==0;
            isRight  = @(coor)  abs(coor(:,1))==xMax;
            isTop    = @(coor)  abs(coor(:,3))==zMax;
            isBottom = @(coor)  abs(coor(:,3))==0;
            isFront  = @(coor)  abs(coor(:,2))==yMax;
            isBack   = @(coor)  abs(coor(:,2))==0;
            isMiddle = @(coor)  abs(coor(:,3))==zMax/2;

            isInSquare = @(coor) (coor(:,1) >= 0.25 & coor(:,1) <= 0.75) & (coor(:,2) >= 0.25 & coor(:,2) <= 0.75);
            
            % Dirichlet

            sDir2.domain    = @(coor) isBottom(coor);
            sDir2.direction = [1,2,3];
            sDir2.value     = 0;
            s.dirichletFun =  DirichletCondition(obj.mesh, sDir2);

            % Neumann
            sPL.domain    = @(coor) isInSquare(coor);
            sPL.direction = 3;
            sPL.value     = -1;
            s.pointloadFun = [];%DistributedLoad(obj.mesh, sPL);

            topFace = obj.mesh.createBoundaryMesh{6};
            bMtop = topFace.mesh;
            originalNodes = topFace.globalConnec(:);
            newNodes      = bMtop.connec(:);
            l2gTop(newNodes(:)) = originalNodes(:);

            [bM,l2g] = bMtop.getBoundarySubmesh(sPL.domain);

%             [bM,l2g] = obj.mesh.getBoundarySubmesh(sPL.domain);

            sAF.fHandle = @(x) [0*x(1,:,:);0*x(3,:,:);sPL.value*ones(size(x(1,:,:)))];
            sAF.ndimf   = 3;
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

            Fext = zeros(bMtop.nnodes,3);
            Fext(l2g,:) = Fext3;

            FextFi = zeros(obj.mesh.nnodes,3);
            FextFi(l2gTop,:) = Fext;
%             aa.fValues = FextFi;
%             aa.mesh = obj.mesh;
%             aa.order = 'P1';
%             p1 = LagrangianFunction(aa);
% %             p1.print('forces3d')
            obj.FextInitial = FextFi; 
            
            s.periodicFun  = [];
            s.mesh         = obj.mesh;

            bc = BoundaryConditions(s);
            obj.boundaryConditions = bc;
        end

    end

end
