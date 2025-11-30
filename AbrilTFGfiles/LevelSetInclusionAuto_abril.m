classdef LevelSetInclusionAuto_abril < handle
    
    properties (Access = public)
        stiffness
        strain
        stress
        dLambda
        bcApplier
        centroids
        nelem
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

        function [obj, u, L, mesh,Kcoarse] = LevelSetInclusionAuto_abril(r, i,nelem,doplot)
            obj.init(r, i,nelem)
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
                  fileName = ['r03_Training' num2str(i)];

                  obj.computeCentroid();
                  CoarsePlotSolution(uFeFun, obj.mesh, obj.bcApplier,fileName, r, obj.centroids);
                  %uFeFun.print(fileName,'Paraview');
                end
             end

             Kcoarse=u.'*obj.stiffness*u;

            
            %% 
            
        end

    end

    methods (Access = private)

        function init(obj, r, i,nelem)
            % close all;
            % clc;
            obj.radius = r;
            obj.nodeDirection = i;
            obj.nelem=nelem;
        end

        function createMesh(obj)
             obj.mesh   = obj.createReferenceMesh();
             %lvSet     = obj.createLevelSetFunction(obj.mesh);
             %uMesh     = obj.computeUnfittedMesh(obj.mesh,lvSet);
             %obj.mesh  = uMesh.createInnerMesh();
             
             obj.boundaryMesh = obj.mesh.createBoundaryMesh();
            [obj.boundaryMeshJoined, obj.localGlobalConnecBd] = obj.mesh.createSingleBoundaryMesh();
        end

        function mesh = createReferenceMesh(obj)
            n       = obj.nelem;
            x1      = linspace(-1,1,n);
            x2      = linspace(-1,1,n);
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

        function [u, L] = doElasticProblemHere(obj)
            s = obj.createElasticProblem();
            obj.createDisplacementFunHere();
            obj.createBCApplyerHere(s);
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
            %obj.computeStrainHere();
            %obj.computeStressHere();

            u=full(u);
            L=full(L);
        end


        function s = createElasticProblem(obj)
            s.mesh     = obj.mesh;
            s.scale    = 'MACRO';
            obj.createMaterial();
            s.material = obj.material;
            s.dim      = '2D';
            s.boundaryConditions = obj.createBoundaryConditions();
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
            E1  = 1;  E2 = E1/1000;
            nu = 1/3;
            x0=0;  y0=0;
            f       = @(x) (sqrt((x(1,:,:)-x0).^2+(x(2,:,:)-y0).^2)<obj.radius)*E2 + ...
                        (sqrt((x(1,:,:)-x0).^2+(x(2,:,:)-y0).^2)>=obj.radius)*E1 ; 
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

            %isForce = @(coor) (abs(coor(:,1) - xMin)   < 1e-10);
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




        function createDisplacementFunHere(obj)
            obj.displacementFun = LagrangianFunction.create(obj.mesh, obj.mesh.ndim, 'P1');
        end

        function createBCApplyerHere(obj, cParams)
            s.mesh = obj.mesh;
            s.boundaryConditions = cParams.boundaryConditions;
            obj.bcApplier = BCApplier(s);
        end

        function createSolverHere(obj, cParams)
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

        %function computeForcesHere(obj, cParams)
        %    n.type         = 'Elastic';
        %    n.scale        = 'MACRO';
        %    n.dim          = obj.getFunDimsHere();
        %    n.BC           = cParams.boundaryConditions;
        %    n.mesh         = obj.mesh;
        %    n.material     = obj.material;
        %    n.globalConnec = obj.mesh.connec;
%
        %    RHSint = RHSIntegrator.create(n);
        %    rhs    = RHSint.compute();
        %    % Perhaps move it inside RHSint?
        %    if strcmp(cParams.solverType,'REDUCED')
        %        R          = RHSint.computeReactions(obj.stiffness);
        %        obj.forces = rhs+R;
        %    else
        %        obj.forces = rhs;
        %    end
        %    
        %end

        function Cg = computeConstraintMatrix(obj)
            obj.dLambda  = LagrangianFunction.create(obj.boundaryMeshJoined, obj.mesh.ndim, 'P1'); 
            test   = LagrangianFunction.create(obj.mesh, obj.mesh.ndim, 'P1');
            f = @(u,v) DP(v,u);
            Cg = IntegrateLHS(f,test,obj.dLambda,obj.mesh,'Boundary',2);   
        end

        function dim = getFunDimsHere(obj)
            d.ndimf     = obj.displacementFun.ndimf;
            d.nnodes    = size(obj.displacementFun.fValues, 1);
            d.ndofs     = d.nnodes*d.ndimf;
            d.nnodeElem = obj.mesh.nnodeElem; % should come from interp..
            d.ndofsElem = d.nnodeElem*d.ndimf;
            dim         = d;
        end


        function [u, L] = computeDisplacementHere(obj,c, rdir)
            K = obj.stiffness;
            nC  = size(c,2);
            Z   = zeros(nC);
            LHS = [K, c; c' Z];

            forces = zeros(obj.displacementFun.nDofs, size(rdir, 2));

            RHS = [forces; rdir];
            sol = LHS\RHS;
            u = sol(1:obj.displacementFun.nDofs,:);
            L = -sol(obj.displacementFun.nDofs+1:end,:); 
%             obj.displacementFun.fValues = u;
%             EIFEMtesting.plotSolution(u,obj.mesh,1,1,0,0)
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

        function rDir = RHSdirichlet(obj) %Son les u_d a sota del vector F
            test   = LagrangianFunction.create(obj.boundaryMeshJoined, obj.mesh.ndim, 'P1');

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

            fB     = {f1x f1y f2x f2y f3x f3y f4x f4y}; 
            nfun = size(fB,2);
            rDir = [];

            for i=1:nfun
                Ud{i}  = AnalyticalFunction.create(fB{i},obj.boundaryMeshJoined);
                f = @(v) DP(v,Ud{i});
                rDire = IntegrateRHS(f,test,obj.boundaryMeshJoined,'Domain',2);
                rDir = [rDir rDire];
            end
        end

        function computeCentroid(obj)
            x0=mean(obj.mesh.coord(:,1));
            y0=mean(obj.mesh.coord(:,2));
            obj.centroids = [x0,y0];
        end

    end
end
