classdef BoundaryConditionsSwan < handle

    properties (Access = private)
        mesh
    end

    methods (Access = public)
        
        function obj = BoundaryConditionsSwan(cParams)
            obj.init(cParams)
        end

         function bc = createBoundaryConditionsTest(obj)
            xMax    = max(obj.mesh.coord(:,1));
            yMax    = max(obj.mesh.coord(:,2));
            isDir   = @(coor)  abs(coor(:,1))==0;
            isForce = @(coor)  (abs(coor(:,1))==xMax & abs(coor(:,2))>=0.499*yMax & abs(coor(:,2))<=0.501*yMax);

            sDir{1}.domain    = @(coor) isDir(coor);
            sDir{1}.direction = [1,2];
            sDir{1}.value     = 0;

            sPL{1}.domain    = @(coor) isForce(coor);
            sPL{1}.direction = 2;
            sPL{1}.value     = -1;

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

         function bc = createBoundaryConditionsTutorialBeam(obj)
            xMax    = max(obj.mesh.coord(:,1));
            yMax    = max(obj.mesh.coord(:,2));
            isDir   = @(coor)  abs(coor(:,1))==0;
            isForce = @(coor)  (abs(coor(:,1))==xMax & abs(coor(:,2))>=0.35*yMax & abs(coor(:,2))<=0.65*yMax);

            sDir{1}.domain    = @(coor) isDir(coor);
            sDir{1}.direction = [1,2];
            sDir{1}.value     = 0;

            sPL{1}.domain    = @(coor) isForce(coor);
            sPL{1}.direction = 2;
            sPL{1}.value     = -1;

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

         function bc = createBoundaryConditionsTutorialBridge(obj)
             xMax    = max(obj.mesh.coord(:,1));
             yMax    = max(obj.mesh.coord(:,2));
             isDir1  = @(coor)  abs(coor(:,2))==0 & abs(coor(:,1))<=0.3;
             isDir2  = @(coor)  abs(coor(:,2))==0 & abs(coor(:,1))>=5.7 & abs(coor(:,1))<=6;
             isForce = @(coor)  abs(coor(:,2))==yMax & abs(coor(:,1))>=2.85 & abs(coor(:,1))<=3.15;

             sDir{1}.domain    = @(coor) isDir1(coor);
             sDir{1}.direction = [2];
             sDir{1}.value     = 0;

             sDir{2}.domain    = @(coor) isDir2(coor);
             sDir{2}.direction = [1,2];
             sDir{2}.value     = 0;

             sPL{1}.domain    = @(coor) isForce(coor);
             sPL{1}.direction = 2;
             sPL{1}.value     = -1;

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

         function bc = createBoundaryConditionsTutorialArch(obj)
             xMax    = max(obj.mesh.coord(:,1));
             yMax    = max(obj.mesh.coord(:,2));
             isDir1  = @(coor)  abs(coor(:,2))==0 & abs(coor(:,1))<=0.2;
             isDir2  = @(coor)  abs(coor(:,2))==0 & abs(coor(:,1))>=1.8 & abs(coor(:,1))<=2;
             isForce = @(coor)  abs(coor(:,2))==0 & abs(coor(:,1))>=0.9 & abs(coor(:,1))<=1.1;

             sDir{1}.domain    = @(coor) isDir1(coor);
             sDir{1}.direction = [1,2];
             sDir{1}.value     = 0;

             sDir{2}.domain    = @(coor) isDir2(coor);
             sDir{2}.direction = [1,2];
             sDir{2}.value     = 0;

             sPL{1}.domain    = @(coor) isForce(coor);
             sPL{1}.direction = 2;
             sPL{1}.value     = -1;

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

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh = cParams.mesh;
        end
    end
end
