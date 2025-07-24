classdef HyperelasticityTestBCs < handle
    
    properties (Access = public)
        boundaryConditions
        dir
    end
    
    properties (Access = private)
        mesh
    end

    methods (Access = public)

        function obj = HyperelasticityTestBCs(type, mesh, perc)
            obj.mesh = mesh;
            switch type
                case 'Traction'
                    obj.createBC_2DTraction(perc);
                case 'Hole'
                    obj.createBC_2DHole();
                case 'HoleDirich'
                    obj.createBC_2DHoleDirich(perc);
                case 'Bending'
                    obj.createBC_2DBending();
                case 'Cube'
                    obj.createBC_3DCube();
                case 'Metamaterial'
                    obj.createBC_2DMetamaterial(perc);
            end
        end

    end

    methods (Access = private)

        function bc = createBC_2DMetamaterial(obj,perc)
            xMax    = max(obj.mesh.coord(:,1));
            yMax    = max(obj.mesh.coord(:,2));
            isLeft   = @(coor)  abs(coor(:,1))==0;
            isRight  = @(coor)  abs(coor(:,1))==xMax;
            isHalf   = @(coor)  abs(abs(coor(:,1)) - xMax/2) <=10e-2;
            isTop    = @(coor)  abs(coor(:,2))==yMax;
            isBottom = @(coor)  abs(coor(:,2))==0;
            isMiddle = @(coor)  abs(abs(coor(:,2))-yMax/2) <= 10e-2;
            
            % 2D N ELEMENTS
            sDir.domain    = @(coor) isLeft(coor);
            sDir.direction = [1,2];
            sDir.value     = 0;
            dir1 =  DirichletCondition(obj.mesh, sDir);

            sDir2.domain    = @(coor) isRight(coor);
            sDir2.direction = [1];
            sDir2.value     = perc*1;
            dir2 =  DirichletCondition(obj.mesh, sDir2);

            sDir3.domain    = @(coor) isRight(coor);
            sDir3.direction = [2];
            sDir3.value     = 0;
            dir3 =  DirichletCondition(obj.mesh, sDir3);




            s.dirichletFun = [dir1, dir2,dir3];
            % 
            % sPL.domain    = @(coor) isRight(coor);
            % sPL.direction = 1;
            % sPL.value     = 0.2;
             s.pointloadFun = [];%DistributedLoad(obj.mesh, sPL);
            % 
            % [bM,l2g] = obj.mesh.getBoundarySubmesh(sPL.domain);
            % 
            % sAF.fHandle = @(x) [sPL.value*ones(size(x(1,:,:)));0*x(2,:,:)];
            % sAF.ndimf   = 2;
            % sAF.mesh    = bM;
            % xFun = AnalyticalFunction(sAF);
            % xFunP1  =xFun.project('P1');
            % 
            % s.mesh = bM;
            % s.type = 'ShapeFunction';
            % s.quadType = 2;
            % rhsI       = RHSintegrator.create(s);
            % test = LagrangianFunction.create(bM,xFun.ndimf,'P1');
            % Fext2 = rhsI.compute(xFunP1,test);   
            % Fext3 = reshape(Fext2,[bM.ndim,bM.nnodes])';
            % 
            % Fext = zeros(obj.mesh.nnodes,2);
            % Fext(l2g,:) = Fext3; 
            
            s.periodicFun  = [];
            s.mesh         = obj.mesh;

            bc = BoundaryConditions(s);
            obj.boundaryConditions = bc;
        end
        function [bc,dir] = createBC_2DTraction(obj,perc)
            xMax    = max(obj.mesh.coord(:,1));
            yMax    = max(obj.mesh.coord(:,2));
            xMin    = min(obj.mesh.coord(:,1));
            yMin    = min(obj.mesh.coord(:,2));

            isLeft   = @(coor)  abs(coor(:,1)-xMin)<= 1e-10;
            isRight  = @(coor)  abs(coor(:,1)-xMax)<= 1e-10;
            isHalf   = @(coor)  abs(abs(coor(:,1)) - xMax/2) <=10e-2;
            isTop    = @(coor)  abs(coor(:,2)-yMax)<= 1e-10;
            isBottom = @(coor)  abs(coor(:,2)-yMin) <= 1e-10;
            isMiddle = @(coor)  abs(abs(coor(:,2))-yMax/2) <= 10e-2;
            
%             % 2D N ELEMENTS
%             sDir.domain    = @(coor) isLeft(coor);
%             sDir.direction = [1,2];
%             sDir.value     = 0;
%             dir1 =  DirichletCondition(obj.mesh, sDir);
% 
%             sDir2.domain    = @(coor) isRight(coor);
%             sDir2.direction = [1];
%             sDir2.value     = perc*6;
%             dir2 =  DirichletCondition(obj.mesh, sDir2);
% 
%             sDir3.domain    = @(coor) isRight(coor);
%             sDir3.direction = [2];
%             sDir3.value     = 0;
%             dir3 =  DirichletCondition(obj.mesh, sDir3);
%             s.dirichletFun = [dir1, dir2,dir3];

 % 2D N ELEMENTS
            sDir{1}.domain    = @(coor) isLeft(coor);
            sDir{1}.direction = [1,2];
            sDir{1}.value     = 0;
%             dir1 =  DirichletCondition(obj.mesh, sDir);

%             sDir{2}.domain    = @(coor) isRight(coor);
%             sDir{2}.direction = [1];
%             sDir{2}.value     = perc*6;
% %             dir2 =  DirichletCondition(obj.mesh, sDir2);

%             sDir{2}.domain    = @(coor) isRight(coor);
%             sDir{2}.direction = [2];
%             sDir{2}.value     = 0;
%             dir3 =  DirichletCondition(obj.mesh, sDir3);

            dirichletFun = [];
            for i = 1:numel(sDir)
                dir = DirichletCondition(obj.mesh, sDir{i});
                dirichletFun = [dirichletFun, dir];
            end
            

            s.dirichletFun = dirichletFun;
%              s.pointloadFun = [];
            % 
            sPL.domain    = @(coor) isRight(coor);
            sPL.direction = [1];
            sPL.value     = 0.15;
             s.pointloadFun = PointLoad(obj.mesh,sPL);%DistributedLoad(obj.mesh, sPL);
            % 
            % [bM,l2g] = obj.mesh.getBoundarySubmesh(sPL.domain);
            % 
            % sAF.fHandle = @(x) [sPL.value*ones(size(x(1,:,:)));0*x(2,:,:)];
            % sAF.ndimf   = 2;
            % sAF.mesh    = bM;
            % xFun = AnalyticalFunction(sAF);
            % xFunP1  =xFun.project('P1');
            % 
            % s.mesh = bM;
            % s.type = 'ShapeFunction';
            % s.quadType = 2;
            % rhsI       = RHSintegrator.create(s);
            % test = LagrangianFunction.create(bM,xFun.ndimf,'P1');
            % Fext2 = rhsI.compute(xFunP1,test);   
            % Fext3 = reshape(Fext2,[bM.ndim,bM.nnodes])';
            % 
            % Fext = zeros(obj.mesh.nnodes,2);
            % Fext(l2g,:) = Fext3; 
            
            s.periodicFun  = [];
            s.mesh         = obj.mesh;

            bc = BoundaryConditions(s);
            obj.boundaryConditions = bc;
            obj.dir =  sDir;
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

%             obj.FextInitial = Fext; 
            
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

%             obj.FextInitial = Fext; 
            
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

%             obj.FextInitial = Fext; 
            
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
%             obj.FextInitial = FextFi; 
            
            s.periodicFun  = [];
            s.mesh         = obj.mesh;

            bc = BoundaryConditions(s);
            obj.boundaryConditions = bc;
        end

    end

end

