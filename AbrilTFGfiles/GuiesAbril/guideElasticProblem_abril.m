classdef guideElasticProblem_abril < handle

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

        function [obj,u,l]=guideElasticProblem_abril(r)
            obj.init(r, 1)
            obj.createMesh();
            [u, L] = obj.doElasticProblemHere();
            u=1;
            l=2;
        end


    end

    methods (Access=private)

        function init (obj,r,i)
            obj.radius=r;
            obj.nodeDirection = i;
        end

        function createMesh(obj)
            bgMesh=obj.createReferenceMesh();  % Crea la background mesh
            lvSet=obj.createlvSetFunction(bgMesh);  % Crea la levelSet per fer el forat
            uMesh=obj.computeUnfittedMesh(bgMesh,lvSet);  %Fa el forat
            obj.mesh = uMesh.createInnerMesh();    % guarda com a mesh el conjunt 
            obj.boundaryMesh = obj.mesh.createBoundaryMesh();   %crea el boundary de la mesh
            [obj.boundaryMeshJoined, obj.localGlobalConnecBd] = obj.mesh.createSingleBoundaryMesh(); %conectivitats
        end


            function bgMesh=createReferenceMesh()
                % Generate coordinates and discretization
                x1 = linspace(-5,5,20);
                x2 = linspace(0,10,20);
            
                % Create the grid
                [xv,yv] = meshgrid(x1,x2);
            
                % Triangulate the mesh to obtain coordinates and connectivities
                [F,V] = mesh2tri(xv,yv,zeros(size(xv)),'x');
                s.coord  = V(:,1:2);
                s.connec = F;
                bgMesh = Mesh.create(s);
            end
    
            function levelSet=createlvSetFunction(obj,bgMesh)
                sLS.type        = 'CircleInclusion';
                sLS.xCoorCenter = 0;    %coordenada y del centre
                sLS.yCoorCenter = 5;  %coordenada y del centre
                sLS.radius      = obj.radius;
                g               = GeometricalFunction(sLS);
                lsFun           = g.computeLevelSetFunction(bgMesh);
                levelSet        = lsFun.fValues;
            end
    
            function uMesh=computeUnfittedMesh(bgMesh,lvSet)
                sUm.backgroundMesh  = bgMesh;          % Defineix la background mesh
                sUm.boundaryMesh    = bgMesh.createBoundaryMesh();   %calcula el boundary mesh corresponent
                uMesh               = UnfittedMesh(sUm);   % background+boundary
                uMesh.compute(lvSet);     % inclusio radi
            end

        function [u,L]=doElasticProblemHere(obj)
            s=obj.createFEMContainer()
            obj.createDisplacementFunHere();
            obj.createBCApplyerHere(s);
            obj.createSolverHere(s)
            obj.computeStiffnessMatrixHere();
            obj.computeForcesHere(s);
             c = obj.computeCmatP1();
             rdir = obj.RHSdir();
            [u, L]  = obj.computeDisplacementHere(c, rdir);
        end


            function s=createFEMContainer(obj)   %Aqui es defineixen les propietats necessaries pel FEM  sobre el solver
                s.mesh     = obj.mesh;  % malla a resoldre
                s.scale    = 'MACRO';   
                obj.createMaterial();
                s.material = obj.material;
                s.dim      = '2D';          % Dimensions del problema
                s.boundaryConditions = obj.createBoundaryConditions();
                s.interpolationType  = 'LINEAR';
                s.solverType         = 'MONOLITHIC';
                s.solverMode         = 'DISP';
                s.solverCase         = 'DIRECT';
            end

            function createMaterial(obj,s)
                [young,poisson] = obj.computeElasticProperties();
                s.type       = 'ISOTROPIC';
                s.ptype      = 'ELASTIC';
                s.ndim       = obj.mesh.ndim;
                s.young      = young;
                s.poisson    = poisson;
                obj.material = Material.create(s);            
            end

            %Aquesta funcio es per crear els 2 materials per tenir el mateix
            %scaling en tots els casos de radi, com en aquest exemple guia
            %nomes es fa per un radi no fa falta i es fara el exemple del
            %tutorial estandard
            function [young,poisson] = computeElasticProperties(obj)
                E1  = 1;
%               E2 = E1/1000;
                nu = 1/3;
%                x0=0;
%                y0=0;
%                 young   = ConstantFunction.create(E,mesh);
%                 poisson = ConstantFunction.create(nu,mesh);
%                f   = @(x) (sqrt((x(1,:,:)-x0).^2+(x(2,:,:)-y0).^2)<obj.radius)*E2 + ...
%                            (sqrt((x(1,:,:)-x0).^2+(x(2,:,:)-y0).^2)>=obj.radius)*E1 ; 
%                                          x(2,:,:).*0 ];
%                young   = AnalyticalFunction.create(f,1,mesh);
                 E       = ConstantFunction.create(E1,obj.mesh);  
                 poisson = ConstantFunction.create(nu,obj.mesh);  
    
            end

            function bc=createBoundaryConditions(obj)
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

            function createDisplacementFunHere(obj)
                obj.uFun = LagrangianFunction.create(obj.mesh, obj.mesh.ndim, 'P1');
            end

            function createBCApplyerHere(obj, cParams)
                s.mesh = obj.mesh;
                s.boundaryConditions = cParams.boundaryConditions;
                obj.bcApplier = BCApplier(s);
            end

            function createSolverHere(obj,cParams)
                sS.type      = cParams.solverCase;
                solver       = Solver.create(sS);
                p.solverType = cParams.solverType;
                p.solverMode = cParams.solverMode;
                p.solver     = solver;
                p.boundaryConditions = cParams.boundaryConditions;
                p.BCApplier          = obj.bcApplier;
                obj.problemSolver    = ProblemSolver(p);
            end

    end
end
