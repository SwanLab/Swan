classdef Tutorial02FEMElasticity < handle

    properties (Access = private)
        mesh
        young
        poisson
        material
        stateProblem
    end

    methods (Access = public)

        function obj = Tutorial02FEMElasticity()
            obj.init();
            obj.createMesh();
            obj.computeElasticProperties();
            obj.createMaterial();
            obj.solveElasticProblem();
        end

    end

    methods (Access = private)

        function init(obj)

        end

        function createMesh(obj)
            obj.mesh = UnitQuadMesh(2,2);
        end

        function computeElasticProperties(obj)
            E1  = 10;
            nu1 = 0.3;
            E   = AnalyticalFunction.create(@(x) E1*ones(size(squeeze(x(1,:,:)))),1,obj.mesh);
            nu  = AnalyticalFunction.create(@(x) nu1*ones(size(squeeze(x(1,:,:)))),1,obj.mesh);
            obj.young   = E;
            obj.poisson = nu;
        end

        function createMaterial(obj)
            s.type    = 'ISOTROPIC';
            s.ptype   = 'ELASTIC';
            s.ndim    = obj.mesh.ndim;
            s.young   = obj.young;
            s.poisson = obj.poisson;
            tensor    = Material.create(s);
            obj.material = tensor;
        end

        function solveElasticProblem(obj)
            s.mesh = obj.mesh;
            s.scale = 'MACRO';
            s.material = obj.material;
            s.dim = '2D';
            s.boundaryConditions = obj.createBoundaryConditions();
            s.solverType = 'MONOLITHIC';
            s.solverMode = 'DISP';
            fem = ElasticProblem(s);
            fem.solve();
            obj.stateProblem = fem;
        end

        function bc = createBoundaryConditions(obj)
%             xMax    = max(obj.mesh.coord(:,1));
% %             yMax    = max(obj.mesh.coord(:,2));
%             isDir   = @(coor)  abs(coor(:,1))==0;
% %             isForce = @(coor)  (abs(coor(:,1))==xMax & abs(coor(:,2))>=0.3*yMax & abs(coor(:,2))<=0.7*yMax);
%             isForce = @(coor)  abs(coor(:,1))==xMax;
% 
%             sDir{1}.domain    = @(coor) isDir(coor);
%             sDir{1}.direction = [1,2];
%             sDir{1}.value     = 0;
% 
%             sPL{1}.domain    = @(coor) isForce(coor);
%             sPL{1}.direction = 2;
%             sPL{1}.value     = -1;
% 
%             dirichletFun = [];
%             for i = 1:numel(sDir)
%                 dir = DirichletCondition(obj.mesh, sDir{i});
%                 dirichletFun = [dirichletFun, dir];
%             end
%             s.dirichletFun = dirichletFun;
% 
%             pointloadFun = [];
%             for i = 1:numel(sPL)
%                 pl = PointLoad(obj.mesh, sPL{i});
%                 pointloadFun = [pointloadFun, pl];
%             end
%             s.pointloadFun = pointloadFun;
% 
%             s.periodicFun  = [];
%             s.mesh = obj.mesh;
%             bc = BoundaryConditions(s);

            xMax    = max(obj.mesh.coord(:,1));
            yMax    = max(obj.mesh.coord(:,2));
            isLeft   = @(coor)  abs(coor(:,1))==0;
            isRight  = @(coor)  abs(coor(:,1))==xMax;
            isTop    = @(coor)  abs(coor(:,2))==yMax;
            isBottom = @(coor)  abs(coor(:,2))==0;
            isMiddle = @(coor)  abs(coor(:,2))==yMax/2;
            
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
            sPL.value     = 0.01;
            s.pointloadFun = PointLoad(obj.mesh, sPL);
            
            s.periodicFun  = [];
            s.mesh         = obj.mesh;

            bc = BoundaryConditions(s);
%             obj.boundaryConditions = bc;

        end

    end

end