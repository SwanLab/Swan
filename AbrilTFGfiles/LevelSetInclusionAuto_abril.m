classdef LevelSetInclusionAuto_abril < handle
    
    properties (Access = public)
        stiffness
        strain
        stress
        dLambda
        bcApplier
    end
    
    properties (Access = private)
        radius
        nodeDirection
        physicalProblem
        
        boundaryConditions
        problemSolver
        forces
        uFun
        strainFun
        
    end

    properties  (Access = protected)
        mesh
        material
        displacementFun
        boundaryMesh
        boundaryMeshJoined
        localGlobalConnecBd
    end


    methods (Access = public)

        function [obj, u, L, mesh,Kcoarse] = LevelSetInclusionAuto_abril(r, i,doplot)
            obj.init(r, i)
            obj.createMesh();
            
            %% New ugly chunk of code warning
            [u, L] = obj.doElasticProblemHere();
            mesh = obj.mesh;
            
             z.mesh      = obj.mesh;
             z.order     = 'P1';

             if doplot==true()
                for i=1:8
                  z.fValues   = reshape(u(:,i),[obj.mesh.ndim,obj.mesh.nnodes])';
                  uFeFun = LagrangianFunction(z);%
                  fileName = ['r05_test' num2str(i)];
                  uFeFun.print(fileName,'Paraview');
                end
             end

             Kcoarse=u.'*obj.stiffness*u;

             %disp(Kcoarse);
            
            %% 
            
        end

    end

    methods (Access = private)

        function init(obj, r, i)
            % close all;
            % clc;
            obj.radius = r;
            obj.nodeDirection = i;
        end

        function createMesh(obj)
            bgMesh   = obj.createReferenceMesh();
             %lvSet    = obj.createLevelSetFunction(bgMesh);
             %uMesh    = obj.computeUnfittedMesh(bgMesh,lvSet);
             %obj.mesh = uMesh.createInnerMesh();
            obj.mesh = bgMesh;
             
            obj.boundaryMesh = obj.mesh.createBoundaryMesh();
            [obj.boundaryMeshJoined, obj.localGlobalConnecBd] = obj.mesh.createSingleBoundaryMesh();
        end

        function mesh = createReferenceMesh(~)
           
             %UnitMesh better
            x1      = linspace(-1,1,20);
            x2      = linspace(-1,1,20);
            [xv,yv] = meshgrid(x1,x2);
            [F,V]   = mesh2tri(xv,yv,zeros(size(xv)),'x');
            s.coord  = V(:,1:2);
            s.connec = F;
            mesh = Mesh.create(s);
        end

        function levelSet = createLevelSetFunction(obj,bgMesh)
            sLS.type        = 'CircleInclusion';
            sLS.xCoorCenter = 0;
            sLS.yCoorCenter = 0;
            sLS.radius      = obj.radius;
            g               = GeometricalFunction(sLS);
            lsFun           = g.computeLevelSetFunction(bgMesh);
            levelSet        = lsFun.fValues;
        end

        function uMesh = computeUnfittedMesh(~,bgMesh,levelSet)
            sUm.backgroundMesh = bgMesh;
            sUm.boundaryMesh   = bgMesh.createBoundaryMesh();
            uMesh              = UnfittedMesh(sUm);
            uMesh.compute(levelSet);
        end

        function s = createElasticProblem(obj)
            s.mesh     = obj.mesh;
            s.scale    = 'MACRO';
            obj.createMaterial();
            s.material = obj.material;
            s.dim      = '2D';
            s.boundaryConditions = obj.createBoundaryConditions();
            %s.interpolationType  = 'LINEAR';
            s.solverType         = 'MONOLITHIC';
            s.solverMode         = 'DISP';
            s.solverCase         = DirectSolver();
        end

        function createMaterial(obj)
            [young,poisson] = obj.computeElasticProperties(obj.mesh);
            s.type       = 'ISOTROPIC';
            s.ptype      = 'ELASTIC';
            s.ndim       = obj.mesh.ndim;
            s.young      = young;
            s.poisson    = poisson;
            obj.material = Material.create(s);
            
        end

        function [young,poisson] = computeElasticProperties(obj,mesh)
            E1  = 1;
            E2 = E1/1000;
            nu = 1/3;
            x0=0;
            y0=0;
%             young   = ConstantFunction.create(E,mesh);
%             poisson = ConstantFunction.create(nu,mesh);
            f   = @(x) (sqrt((x(1,:,:)-x0).^2+(x(2,:,:)-y0).^2)<obj.radius)*E2 + ...
                        (sqrt((x(1,:,:)-x0).^2+(x(2,:,:)-y0).^2)>=obj.radius)*E1 ; 
%                                      x(2,:,:).*0 ];
            young   = AnalyticalFunction.create(f,mesh);
            poisson = ConstantFunction.create(nu,mesh);            
        end

        function bc = createBoundaryConditions(obj)

            v                    = zeros(8,1);
            v(obj.nodeDirection) = 1;
            nRes                 = [1 1 2 2 3 3 4 4]*v;
            assignMatrix         = [2 1 0 0 0 0 0 0
                                    0 0 2 1 0 0 0 0
                                    0 0 0 0 2 1 0 0
                                    0 0 0 0 0 0 2 1
                                    1 2 1 2 1 2 1 2];

            vSimp       = assignMatrix*v;
            dirs        = cell(5,1);
            [dirs{:,1}] = deal([1, 2]);
            [dirs{nRes}]  = deal(vSimp(nRes));
            [dirs{end}]   = deal(vSimp(end));

            xMax    = max(obj.mesh.coord(:,1));
            yMax    = max(obj.mesh.coord(:,2));
            xMin    = min(obj.mesh.coord(:,1));
            yMin    = min(obj.mesh.coord(:,2));
            tol     = 1e-10;

            corner1 = @(coor) abs(coor(:,1)-xMin) <= tol & abs(coor(:,2)-yMin)<= tol;
            corner2 = @(coor) abs(coor(:,1)-xMax) <= tol & abs(coor(:,2)-yMin)<= tol;
            corner3 = @(coor) abs(coor(:,1)-xMin) <= tol & abs(coor(:,2)-yMax)<= tol;
            corner4 = @(coor) abs(coor(:,1)-xMax) <= tol & abs(coor(:,2)-yMax)<= tol;

            cornerVec = {corner1; corner2; corner3; corner4};

            isForce = @(coor) (abs(coor(:,1) - xMin)   < 1e-10);


            sDir{1}.domain    = @(coor) cornerVec{1}(coor);
            sDir{1}.direction = cell2mat(dirs(1));
            sDir{1}.value     = 0;

            
            sDir{2}.domain    = @(coor) cornerVec{2}(coor);
            sDir{2}.direction = cell2mat(dirs(2));
            sDir{2}.value     = 0;

            sDir{3}.domain    = @(coor) cornerVec{3}(coor);
            sDir{3}.direction = cell2mat(dirs(3));
            sDir{3}.value     = 0;

            sDir{4}.domain    = @(coor) cornerVec{4}(coor);
            sDir{4}.direction = cell2mat(dirs(4));
            sDir{4}.value     = 0;

            sDir{5}.domain    = @(coor) cornerVec{nRes}(coor);
            sDir{5}.direction = cell2mat(dirs(end));
            sDir{5}.value     = 1;

            % sPL{1}.domain    = @(coor) isForce(coor);
            % sPL{1}.direction = 2;
            % sPL{1}.value     = 0;
            sPL = {};

            dirichletFun = [];
            for i = 1:numel(sDir)
                dir = DirichletCondition(obj.mesh, sDir{i});
                dirichletFun = cat(2, dirichletFun, dir);
            end
            s.dirichletFun = dirichletFun;

            pointloadFun = [];
            for i = 1:numel(sPL)
                pl = PointLoad(obj.mesh, sPL{i});
                pointloadFun = cat(2, pointloadFun, pl);
            end
            s.pointloadFun = pointloadFun;

            s.periodicFun  = [];
            s.mesh = obj.mesh;
            bc = BoundaryConditions(s);
        end

%         function bc = createBoundaryConditions(obj)
%             Lx = max(obj.mesh.coord(:,1)) - min(obj.mesh.coord(:,1));
%             Ly = max(obj.mesh.coord(:,2)) - min(obj.mesh.coord(:,2));
%                 f1 = @(x) 1/(4*Lx/2*Ly/2)*(Lx/2 - x(1,:,:)+0.5).*(Ly/2 - x(2,:,:)+0.5);
%                 f2 = @(x) 1/(4*Lx/2*Ly/2)*(Lx/2 + x(1,:,:)+0.5).*(Ly/2 - x(2,:,:)+0.5);
%                 f3 = @(x) 1/(4*Lx/2*Ly/2)*(Lx/2 + x(1,:,:)+0.5).*(Ly/2 + x(2,:,:)+0.5);
%                 f4 = @(x) 1/(4*Lx/2*Ly/2)*(Lx/2 - x(1,:,:)+0.5).*(Ly/2 + x(2,:,:)+0.5);
%                 ndimf = 1;
% 
%                 N11 = AnalyticalFunction.create(f1,ndimf,obj.boundaryMesh{1}.mesh);
%                 N13 = AnalyticalFunction.create(f1,ndimf,obj.boundaryMesh{3}.mesh);
%                 N23 = AnalyticalFunction.create(f1,ndimf,obj.boundaryMesh{3}.mesh);
%                 N22 = AnalyticalFunction.create(f1,ndimf,obj.boundaryMesh{2}.mesh);
%                 N23 = AnalyticalFunction.create(f1,ndimf,obj.boundaryMesh{3}.mesh);
%                 N22 = AnalyticalFunction.create(f1,ndimf,obj.boundaryMesh{2}.mesh);
% 
%                 coord = [obj.boundaryMesh{1}.mesh.coord]
% 
% 
%         end

        function [u, L] = doElasticProblemHere(obj)
            s = obj.createElasticProblem();
            obj.createDisplacementFunHere();
            %obj.createBCApplyerHere(s);
            obj.createSolverHere(s)
            obj.computeStiffnessMatrix();
            %obj.computeForcesHere(s);
            c = obj.computeConstraintMatrix();
            rdir = obj.RHSdirichlet();
            [u, L]  = obj.computeDisplacementHere(c, rdir);

            if isa(obj.dLambda, "LagrangianFunction")
                l2g_dof = ((obj.localGlobalConnecBd*obj.displacementFun.ndimf)' - ((obj.displacementFun.ndimf-1):-1:0))';
                l2g_dof = l2g_dof(:);
                uB = u(l2g_dof, :);
                L = uB'*L;
            end
            % obj.plotSolution(u,L);
            %obj.computeStrainHere();
            %obj.computeStressHere();

            u=full(u);
            L=full(L);
        end

        function createDisplacementFunHere(obj)
            obj.displacementFun = LagrangianFunction.create(obj.mesh, obj.mesh.ndim, 'P1');
        end

        function createBCApplyerHere(obj, cParams)
            s.mesh = obj.mesh;
            s.boundaryConditions = cParams.boundaryConditions;
            obj.bcApplier = BCApplier(s);
        end

        function createSolverHere(obj, cParams)
           % sS.type      = cParams.solverCase;
            p.solverType = cParams.solverType;
            p.solverMode = cParams.solverMode;
            p.solver     = cParams.solverCase;

            p.boundaryConditions = cParams.boundaryConditions;
            p.BCApplier          = obj.bcApplier;
            obj.problemSolver    = ProblemSolver(p);
        end

        function computeStiffnessMatrix(obj)
            C     = obj.material;
            f = @(u,v) DDP(SymGrad(v),DDP(C,SymGrad(u)));
            obj.stiffness = IntegrateLHS(f,obj.displacementFun,obj.displacementFun,obj.mesh,'Domain',2);
        end

        function computeForcesHere(obj, cParams)
            n.type         = 'Elastic';
            n.scale        = 'MACRO';
            n.dim          = obj.getFunDimsHere();
            n.BC           = cParams.boundaryConditions;
            n.mesh         = obj.mesh;
            n.material     = obj.material;
            n.globalConnec = obj.mesh.connec;

            RHSint = RHSIntegrator.create(n);
            rhs    = RHSint.compute();
            % Perhaps move it inside RHSint?
            if strcmp(cParams.solverType,'REDUCED')
                R          = RHSint.computeReactions(obj.stiffness);
                obj.forces = rhs+R;
            else
                obj.forces = rhs;
            end
            
        end

        function Cg = computeConstraintMatrix(obj)
            s.quadType = 2;
            s.boundaryMeshJoined    = obj.boundaryMeshJoined;
            s.localGlobalConnecBd   = obj.localGlobalConnecBd;
            s.nnodes                 = obj.mesh.nnodes;

            lhs = LHSintegrator_MassBoundary_albert(s);
            test   = LagrangianFunction.create(obj.boundaryMeshJoined, obj.mesh.ndim, 'P1'); % !!
            obj.dLambda  = LagrangianFunction.create(obj.boundaryMeshJoined, obj.mesh.ndim, 'P1');
            Cg = lhs.compute(obj.dLambda,test); 

            test   = LagrangianFunction.create(obj.mesh, obj.mesh.ndim, 'P1');
            f = @(u,v) DP(v,u);
            Cg2 = IntegrateLHS(f,test,obj.dLambda,obj.mesh,'Boundary',2);   

        end

        function Cg = computeCmatP2(obj)
            s.quadType = 2;
            s.boundaryMeshJoined    = obj.boundaryMeshJoined;
            s.localGlobalConnecBd   = obj.localGlobalConnecBd;
            s.nnodes                 = obj.mesh.nnodes;

            % lhs = LHSintegrator_ShapeFunction_fun(s);
            lhs = LHSintegrator_MassBoundary_albert(s);
            test   = LagrangianFunction.create(obj.boundaryMeshJoined, obj.mesh.ndim, 'P1'); % !!
             obj.dLambda  = LagrangianFunction.create(obj.boundaryMeshJoined, obj.mesh.ndim, 'P1');
             Cg = lhs.compute(obj.dLambda,test);      
        end

        function Cg = computeCmat(obj)
            s.quadType = 2;
            s.boundaryMeshJoined    = obj.boundaryMeshJoined;
            s.localGlobalConnecBd   = obj.localGlobalConnecBd;
            s.nnodes                 = obj.mesh.nnodes;
            s.mesh = obj.boundaryMeshJoined;

            lhs = LHSintegrator_ShapeFunction_fun(s);
            % lhs = LHSintegrator_MassBoundary_albert(s);
            test   = LagrangianFunction.create(obj.boundaryMeshJoined, obj.mesh.ndim, 'P1'); % !!
            ndimf  = 2; 
            Lx     = max(obj.mesh.coord(:,1)) - min(obj.mesh.coord(:,1));
            Ly     = max(obj.mesh.coord(:,2)) - min(obj.mesh.coord(:,2));
%             f1 = @(x) [1/(4*Lx/2*Ly/2)*(1-x(1,:,:)).*(1-x(2,:,:));...
%                     1/(4*Lx/2*Ly/2)*(1-x(1,:,:)).*(1-x(2,:,:))  ];
%             f2 = @(x) [1/(4*Lx/2*Ly/2)*(1+x(1,:,:)).*(1-x(2,:,:));...
%                     1/(4*Lx/2*Ly/2)*(1+x(1,:,:)).*(1-x(2,:,:))  ];
%             f3 = @(x) [1/(4*Lx/2*Ly/2)*(1+x(1,:,:)).*(1+x(2,:,:));...
%                     1/(4*Lx/2*Ly/2)*(1+x(1,:,:)).*(1+x(2,:,:))  ];
%             f4 = @(x) [1/(4*Lx/2*Ly/2)*(1-x(1,:,:)).*(1+x(2,:,:));...
%                     1/(4*Lx/2*Ly/2)*(1-x(1,:,:)).*(1+x(2,:,:))  ];
            f1 = @(x) [1/(4)*(1-x(1,:,:)).*(1-x(2,:,:));...
                    1/(4)*(1-x(1,:,:)).*(1-x(2,:,:))  ];
            f2 = @(x) [1/(4)*(1+x(1,:,:)).*(1-x(2,:,:));...
                    1/(4)*(1+x(1,:,:)).*(1-x(2,:,:))  ];
            f3 = @(x) [1/(4)*(1+x(1,:,:)).*(1+x(2,:,:));...
                    1/(4)*(1+x(1,:,:)).*(1+x(2,:,:))  ];
            f4 = @(x) [1/(4)*(1-x(1,:,:)).*(1+x(2,:,:));...
                    1/(4)*(1-x(1,:,:)).*(1+x(2,:,:))  ];
            f     = {f1 f2 f3 f4}; %
            nfun = size(f,2);
            Cg = [];
            for i=1:nfun
                obj.dLambda{i}  = AnalyticalFunction.create(f{i},ndimf,obj.boundaryMeshJoined);
                    
                %% Project to P1
%                 obj.dLambda{i} = obj.dLambda{i}.project('P1');

                Ce = lhs.compute(obj.dLambda{i},test);
                [iLoc,jLoc,vals] = find(Ce);
    
%                 l2g_dof = ((obj.localGlobalConnecBd*test.ndimf)' - ((test.ndimf-1):-1:0))';
%                 l2g_dof = l2g_dof(:);
%                 jGlob = l2g_dof(jLoc);
%                 Cg = [Cg sparse(iLoc,jGlob,vals, obj.displacementFun.nDofs, dLambda.nDofs)];

                l2g_dof = ((obj.localGlobalConnecBd*test.ndimf)' - ((test.ndimf-1):-1:0))';
                l2g_dof = l2g_dof(:);
                iGlob = l2g_dof(iLoc);
                Cg = [Cg sparse(iGlob,jLoc,vals, obj.displacementFun.nDofs, obj.dLambda{i}.ndimf)];
            end

        end

        function Cg = computeCmatEdgeCorner(obj,data)
            s.quadType = 2;
            s.mesh     = obj.boundaryMeshJoined;
            lhs = LHSintegrator_ShapeFunction_fun(s);
            test   = LagrangianFunction.create(obj.boundaryMeshJoined, obj.mesh.ndim, 'P1'); % !!
            ndimf  = 2;
            Lx     = max(obj.mesh.coord(:,1)) - min(obj.mesh.coord(:,1));
            Ly     = max(obj.mesh.coord(:,2)) - min(obj.mesh.coord(:,2));
%             f1 = @(x) [1/(4*Lx/2*Ly/2)*(1-x(1,:,:)).*(1-x(2,:,:));...
%                     1/(4*Lx/2*Ly/2)*(1-x(1,:,:)).*(1-x(2,:,:))  ];
%             f2 = @(x) [1/(4*Lx/2*Ly/2)*(1+x(1,:,:)).*(1-x(2,:,:));...
%                     1/(4*Lx/2*Ly/2)*(1+x(1,:,:)).*(1-x(2,:,:))  ];
%             f3 = @(x) [1/(4*Lx/2*Ly/2)*(1+x(1,:,:)).*(1+x(2,:,:));...
%                     1/(4*Lx/2*Ly/2)*(1+x(1,:,:)).*(1+x(2,:,:))  ];
%             f4 = @(x) [1/(4*Lx/2*Ly/2)*(1-x(1,:,:)).*(1+x(2,:,:));...
%                     1/(4*Lx/2*Ly/2)*(1-x(1,:,:)).*(1+x(2,:,:))  ];
            f1 = @(x) [1/(4)*(1-x(1,:,:)).*(1-x(2,:,:));...
                    1/(4)*(1-x(1,:,:)).*(1-x(2,:,:))  ];
            f2 = @(x) [1/(4)*(1+x(1,:,:)).*(1-x(2,:,:));...
                    1/(4)*(1+x(1,:,:)).*(1-x(2,:,:))  ];
            f3 = @(x) [1/(4)*(1+x(1,:,:)).*(1+x(2,:,:));...
                    1/(4)*(1+x(1,:,:)).*(1+x(2,:,:))  ];
            f4 = @(x) [1/(4)*(1-x(1,:,:)).*(1+x(2,:,:));...
                    1/(4)*(1-x(1,:,:)).*(1+x(2,:,:))  ];
            f     = {f1 f2 f3 f4}; %
            nfun = size(f,2);
            Cg = [];
            for i=1:nfun
                obj.dLambda{i}  = AnalyticalFunction.create(f{i},ndimf,obj.boundaryMeshJoined);
                    
                %% Project to P1
%                  obj.dLambda{i} = obj.dLambda{i}.project('P1');

                Ce = lhs.compute(obj.dLambda{i},test);
                [iLoc,jLoc,vals] = find(Ce);
    
%                 l2g_dof = ((obj.localGlobalConnecBd*test.ndimf)' - ((test.ndimf-1):-1:0))';
%                 l2g_dof = l2g_dof(:);
%                 jGlob = l2g_dof(jLoc);
%                 Cg = [Cg sparse(iLoc,jGlob,vals, obj.displacementFun.nDofs, dLambda.nDofs)];

                l2g_dof = ((obj.localGlobalConnecBd*test.ndimf)' - ((test.ndimf-1):-1:0))';
                l2g_dof = l2g_dof(:);
                iGlob = l2g_dof(iLoc);
                Cg = [Cg sparse(iGlob,jLoc,vals, obj.displacementFun.nDofs, obj.dLambda{i}.ndimf)];
            end
            cornerDof = data.boundaryConditions.dirichlet_dofs;
            Cg(cornerDof,:) = 0;
            Cg2 = zeros(size(Cg,1),size(cornerDof,1));
            for i=1:length(cornerDof)
               Cg2(cornerDof(i),i)=1; 
            end
            Cg = [Cg,Cg2];
        end

            function Cg = computeCmatRB(obj,data)
            s.quadType = 2;
            s.mesh     = obj.boundaryMeshJoined;
            lhs = LHSintegrator_ShapeFunction_fun(s);
            test   = LagrangianFunction.create(obj.boundaryMeshJoined, obj.mesh.ndim, 'P1'); % !!
            ndimf  = 2;
            minx   = min(obj.mesh.coord(:,1));
            miny   = min(obj.mesh.coord(:,2));
            maxx   = max(obj.mesh.coord(:,1));
            maxy   = max(obj.mesh.coord(:,2));
            x0     = (maxx+minx)/2;
            y0     = (maxy+miny)/2;
            Lx     = max(obj.mesh.coord(:,1)) - min(obj.mesh.coord(:,1));
            Ly     = max(obj.mesh.coord(:,2)) - min(obj.mesh.coord(:,2));
            theta  = 1*pi/180;
            cond   = [minx, maxx, miny, maxy];
            nrb=3;
            k=1;
            dof  = [1,1,2,2];
            for i=1:length(obj.boundaryMesh)
                x0 = (max(obj.boundaryMesh{i}.mesh.coord(:,1))+ min(obj.boundaryMesh{i}.mesh.coord(:,1)))/2;
                y0 = (max(obj.boundaryMesh{i}.mesh.coord(:,2))+ min(obj.boundaryMesh{i}.mesh.coord(:,2)))/2;
                    f{k}   = @(x) [(x(dof(i),:,:)==cond(i)).*(x(1,:,:)./x(1,:,:)) + (x(dof(i),:,:)~=cond(i)).*(x(1,:,:)*0);...
                                     x(2,:,:).*0 ];
                    f{k+1} = @(x) [x(dof(i),:,:).*0 ;...
                                  (x(dof(i),:,:)==cond(i)).*(x(1,:,:)./x(1,:,:)) + (x(dof(i),:,:)~=cond(i)).*(x(1,:,:)*0) ];
                    f{k+2} = @(x) [(x(dof(i),:,:)==cond(i)).*((x(1,:,:)-x0).*cos(theta) - (x(2,:,:)-y0).*sin(theta)) + (x(dof(i),:,:)~=cond(i)).*(x(1,:,:)*0) ;...
                                   (x(dof(i),:,:)==cond(i)).*((x(1,:,:)-x0).*sin(theta) + (x(2,:,:)-y0).*cos(theta))*-1 + (x(dof(i),:,:)~=cond(i)).*(x(1,:,:)*0) ];
%                     f{k+2} = @(x) [(x(dof(i),:,:)==cond(i)).*((x(2,:,:)-y0)) + (x(dof(i),:,:)~=cond(i)).*(x(1,:,:)*0) ;...
%                                    (x(dof(i),:,:)==cond(i)).*(-(x(1,:,:)-x0)) + (x(dof(i),:,:)~=cond(i)).*(x(1,:,:)*0) ];
                    k=k+3;
            end            
%             f1 = @(x) [(x(2,:,:)==miny).*(x(1,:,:)./x(1,:,:)) + (x(2,:,:)~=miny).*(x(1,:,:)*0);...
%                         x(2,:,:).*0 ];
%             f2 = @(x) [x(2,:,:).*0 ;...
%                         (x(2,:,:)==miny).*(x(1,:,:)./x(1,:,:)) + (x(2,:,:)~=miny).*(x(1,:,:)*0) ];
%              f3 = @(x) [(x(2,:,:)==miny).*((x(1,:,:)-x0).*cos(theta) - (x(2,:,:)-(-1)).*sin(theta)) + (x(2,:,:)~=miny).*(x(1,:,:)*0) ;...
%                         (x(2,:,:)==miny).*((x(1,:,:)-x0).*sin(theta) + (x(2,:,:)-(-1)).*cos(theta)) + (x(2,:,:)~=miny).*(x(1,:,:)*0) ];
% %             f2 = @(x) [1/(4)*(1+x(1,:,:)).*(1-x(2,:,:));...
% %                     1/(4)*(1+x(1,:,:)).*(1-x(2,:,:))  ];
% %             f3 = @(x) [1/(4)*(1+x(1,:,:)).*(1+x(2,:,:));...
% %                     1/(4)*(1+x(1,:,:)).*(1+x(2,:,:))  ];
%             f4 = @(x) [1/(4)*(1-x(1,:,:)).*(1+x(2,:,:));...
%                     1/(4)*(1-x(1,:,:)).*(1+x(2,:,:))  ];
%             f     = {f1 f2 f3 f4}; %
            nfun = size(f,2);
            Cg = [];
            for i=1:nfun
                obj.dLambda{i}  = AnalyticalFunction.create(f{i},ndimf,obj.boundaryMeshJoined);
%                 dLambda = dLambda.project('P1');
%                 obj.dLambda{i}.plot
                Ce = lhs.compute(obj.dLambda{i},test);
                Ce = sum(Ce,2);
                [iLoc,jLoc,vals] = find(Ce);
    
%                 l2g_dof = ((obj.localGlobalConnecBd*test.ndimf)' - ((test.ndimf-1):-1:0))';
%                 l2g_dof = l2g_dof(:);
%                 jGlob = l2g_dof(jLoc);
%                 Cg = [Cg sparse(iLoc,jGlob,vals, obj.displacementFun.nDofs, dLambda.nDofs)];

                l2g_dof = ((obj.localGlobalConnecBd*test.ndimf)' - ((test.ndimf-1):-1:0))';
                l2g_dof = l2g_dof(:);
                iGlob = l2g_dof(iLoc);
                Cg = [Cg sparse(iGlob,jLoc,vals, obj.displacementFun.nDofs, 1)];
            end
%             cornerDof = data.boundaryConditions.dirichlet_dofs;
%             Cg(cornerDof,:) = 0;
%             Cg2 = zeros(size(Cg,1),size(cornerDof,1));
%             for i=1:length(cornerDof)
%                Cg2(cornerDof(i),i)=1; 
%             end
%             Cg = [Cg,Cg2];
        end

     

        function dim = getFunDimsHere(obj)
            d.ndimf     = obj.displacementFun.ndimf;
            d.nnodes    = size(obj.displacementFun.fValues, 1);
            d.ndofs     = d.nnodes*d.ndimf;
            d.nnodeElem = obj.mesh.nnodeElem; % should come from interp..
            d.ndofsElem = d.nnodeElem*d.ndimf;
            dim         = d;
        end

%         function [u, L] = computeDisplacementHere(obj)
%             o.stiffness = obj.stiffness;
%             o.forces    = obj.forces;
%             [u,L]       = obj.problemSolver.solve(o);
%             z.mesh      = obj.mesh;
%             z.fValues   = reshape(u,[obj.mesh.ndim,obj.mesh.nnodes])';
%             z.order     = 'P1';
%             obj.uFun    = LagrangianFunction(z);
%             uSplit      = reshape(u,[obj.mesh.ndim,obj.mesh.nnodes])';
%             obj.displacementFun.fValues = uSplit;
%         end

        function [u, L] = computeDisplacementHere(obj,c, rdir)
            K = obj.stiffness;
            nC  = size(c,2);
            Z   = zeros(nC);
            LHS = [K, c; c' Z];
            % ud = eye(nC);
% %             ud(7) = 1;
%             ud(1) = 1;
            forces = zeros(obj.displacementFun.nDofs, size(rdir, 2));

            RHS = [forces; rdir];
            sol = LHS\RHS;
            u = sol(1:obj.displacementFun.nDofs,:);
            L = -sol(obj.displacementFun.nDofs+1:end,:); 
%             obj.displacementFun.fValues = u;
%             EIFEMtesting.plotSolution(u,obj.mesh,1,1,0,0)
        end

        function plotSolution(obj,u,L)
%             obj.displacementFun.fValues = u;
%             EIFEMtesting.plotSolution(u,obj.mesh,1,1,0,0)
             L = reshape(L,2,[])';
            nbasis = size(obj.dLambda,2);
            for i = 1:nbasis
                lambda{i} = obj.dLambda{i}.*L(i,:)';
                lambdaP1 = lambda{i}.project('P1');
                EIFEMtesting.plotSolution(lambdaP1.fValues,obj.boundaryMeshJoined,2,1,i,0)
            end
        end

        function computeStrainHere(obj)
            quad = Quadrature.create(obj.mesh, 2);
            xV   = quad.posgp;
            obj.strainFun  = SymGrad(obj.displacementFun);
%             strFun       = strFun.obtainVoigtFormat();
            obj.strain     = obj.strainFun.evaluate(xV);
        end

        function computeStressHere(obj)
            quad            = Quadrature.create(obj.mesh, 2);
            xV              = quad.posgp;
            stressFun       = DDP(obj.material, obj.strainFun);
            stressFun.ndimf = obj.strainFun.ndimf;
            obj.stress      = stressFun.evaluate(xV);

        end


        function rDir = RHSdirichlet(obj)
            test   = LagrangianFunction.create(obj.boundaryMeshJoined, obj.mesh.ndim, 'P1');
           % s.mesh = obj.boundaryMeshJoined;
           % s.quadType = 2;
        %    rhs = RHSIntegratorShapeFunction(s);

            f1x = @(x) [1/(4)*(1-x(1,:,:)).*(1-x(2,:,:));...
                    0*x(2,:,:)  ];
            f2x = @(x) [1/(4)*(1+x(1,:,:)).*(1-x(2,:,:));...
                    0*x(2,:,:)  ];
            f3x = @(x) [1/(4)*(1+x(1,:,:)).*(1+x(2,:,:));...
                    0*x(2,:,:)  ];
            f4x = @(x) [1/(4)*(1-x(1,:,:)).*(1+x(2,:,:));...
                    0*x(2,:,:)  ];

            f1y = @(x) [0*x(1,:,:);...
                    1/(4)*(1-x(1,:,:)).*(1-x(2,:,:))  ];
            f2y = @(x) [0*x(1,:,:);...
                    1/(4)*(1+x(1,:,:)).*(1-x(2,:,:))  ];
            f3y = @(x) [0*x(1,:,:);...
                    1/(4)*(1+x(1,:,:)).*(1+x(2,:,:))  ];
            f4y = @(x) [0*x(1,:,:);...
                    1/(4)*(1-x(1,:,:)).*(1+x(2,:,:))  ];


            fB     = {f1x f1y f2x f2y f3x f3y f4x f4y}; %
            nfun = size(fB,2);
            rDir = [];
            for i=1:nfun
                Ud{i}  = AnalyticalFunction.create(fB{i},obj.boundaryMeshJoined);
                    
                % %% Project to P1
                % obj.dLambda{i} = obj.dLambda{i}.project('P1');


          %      rDire = rhs.compute(Ud{i},test);
                 f = @(v) DP(v,Ud{i});
                 rDire = IntegrateRHS(f,test,obj.boundaryMeshJoined,'Domain',2);
%                 [iLoc,jLoc,vals] = find(Ce);
% 
% %                 l2g_dof = ((obj.localGlobalConnecBd*test.ndimf)' - ((test.ndimf-1):-1:0))';
% %                 l2g_dof = l2g_dof(:);
% %                 jGlob = l2g_dof(jLoc);
% %                 Cg = [Cg sparse(iLoc,jGlob,vals, obj.displacementFun.nDofs, dLambda.nDofs)];
% 
%                 l2g_dof = ((obj.localGlobalConnecBd*test.ndimf)' - ((test.ndimf-1):-1:0))';
%                 l2g_dof = l2g_dof(:);
%                 iGlob = l2g_dof(iLoc);
                rDir = [rDir rDire];
            end

           


        end

         function rDir = RHSweak(obj, c)
                rDir = eye(size(c,2));

         end

    end
end
