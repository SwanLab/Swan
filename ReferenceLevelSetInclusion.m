classdef ReferenceLevelSetInclusion < handle

    properties (Access = private)
        mesh
        physicalProblem        
    end

    methods (Access = public)

        function obj = ReferenceLevelSetInclusion()
            obj.init()
            obj.createMesh();
            obj.createElasticProblem();
            obj.physicalProblem.solve();
        end

    end

    methods (Access = private)

        function init(obj)
            close all;
        end

        function createMesh(obj)
            bgMesh = obj.createReferenceMesh();
            lvSet  = obj.createLevelSetFunction(bgMesh);
            uMesh  = obj.computeUnfittedMesh(bgMesh,lvSet);
            obj.mesh = uMesh.createInnerMesh();
        end

        function mesh = createReferenceMesh(obj)
             %UnitMesh better
            x1      = linspace(0,1,50);
            x2      = linspace(0,1,50);
            [xv,yv] = meshgrid(x1,x2);
            [F,V]   = mesh2tri(xv,yv,zeros(size(xv)),'x');
            s.coord  = V(:,1:2);
            s.connec = F;
            mesh = Mesh.create(s);
        end

        function levelSet = createLevelSetFunction(obj,bgMesh)
            sLS.type        = 'CircleInclusion';
            sLS.xCoorCenter = 0.5;
            sLS.yCoorCenter = 0.5;
            sLS.radius      = 0.2;
            g               = GeometricalFunction(sLS);
            lsFun           = g.computeLevelSetFunction(bgMesh);
            levelSet        = lsFun.fValues;
        end

        function uMesh = computeUnfittedMesh(obj,bgMesh,levelSet)
            sUm.backgroundMesh = bgMesh;
            sUm.boundaryMesh   = bgMesh.createBoundaryMesh();
            uMesh = UnfittedMesh(sUm);
            uMesh.compute(levelSet);
        end

        function createElasticProblem(obj)
            s.mesh = obj.mesh;
            s.scale = 'MACRO';
            s.material = obj.createMaterial();
            s.dim = '2D';
            s.boundaryConditions = obj.createBoundaryConditions();
            s.interpolationType = 'LINEAR';
            s.solverType = 'REDUCED';
            s.solverMode = 'DISP';
            s.solverCase = 'DIRECT';
            fem = ElasticProblem(s);
            obj.physicalProblem = fem;
        end

        function material = createMaterial(obj)
            [young,poisson] = obj.computeElasticProperties(obj.mesh);
            s.type    = 'ISOTROPIC';
            s.ptype   = 'ELASTIC';
            s.ndim    = obj.mesh.ndim;
            s.young   = young;
            s.poisson = poisson;
            tensor    = Material.create(s);
            material  = tensor;
        end

        function [young,poisson] = computeElasticProperties(obj,mesh)
            E  = 1;
            nu = 1/3;
            young   = ConstantFunction.create(E,mesh);
            poisson = ConstantFunction.create(nu,mesh);            
        end

        function bc = createBoundaryConditions(obj)
            xMax    = max(obj.mesh.coord(:,1));
            yMax    = max(obj.mesh.coord(:,2));
            xMin    = min(obj.mesh.coord(:,1));
            yMin    = min(obj.mesh.coord(:,2));
            tol     = 1e-10;
            isLeft  = @(coor)  abs(coor(:,1)-xMin) <= tol;
            isRight = @(coor)  abs(coor(:,1)-xMax) <= tol;
%             isForce = @(coor)  (abs(coor(:,1))==xMax & abs(coor(:,2))>=0.4*yMax & abs(coor(:,2))<=0.6*yMax);

            sDir{1}.domain    = @(coor) isLeft(coor);
            sDir{1}.direction = [1,2];
            sDir{1}.value     = 0;

            sDir{2}.domain    = @(coor) isRight(coor);
            sDir{2}.direction = 1;
            sDir{2}.value     = 1;

%             sPL{1}.domain    = @(coor) isForce(coor);
%             sPL{1}.direction = 2;
%             sPL{1}.value     = -1;
            sPL = {};

            dirichletFun = [];
            for i = 1:numel(sDir)
                dir = DirichletCondition(obj.mesh, sDir{i});
                dirichletFun = [dirichletFun, dir];
            end
            s.dirichletFun = dirichletFun;

            pointloadFun = [];
            for i = 1:numel(sPL)
                pl = PointLoad(obj.mesh, sPL{i});
                pointloadFun = [pointloadFun, pl];
            end
            s.pointloadFun = pointloadFun;

            s.periodicFun  = [];
            s.mesh = obj.mesh;
            bc = BoundaryConditions(s);
        end
    end
end