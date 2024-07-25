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
            obj.mesh = Mesh.createFromGiD('hole_mesh.m');
        end

        function computeElasticProperties(obj)
            E1  = 2.909090909090909e+03;
            nu1 = 0.454545454545455;
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
            s.solverType = 'REDUCED';
            s.solverMode = 'DISP';
            s.solverCase = 'DIRECT';
            fem = ElasticProblem(s);
            fem.solve();
            obj.stateProblem = fem;
        end

        function bc = createBoundaryConditions(obj)
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
%             s.dirichletFun = [dir2];

            sDir3.domain    = @(coor) isTop(coor);
            sDir3.direction = [1];
            sDir3.value     = 0;
            dir3 =  DirichletCondition(obj.mesh, sDir3);

            sDir4.domain    = @(coor) isTop(coor);
            sDir4.direction = [2];
            sDir4.value     = 0.5;
            dir4 =  DirichletCondition(obj.mesh, sDir4);

            s.dirichletFun = [dir2,dir3,dir4];
            s.pointloadFun = [];
% 
%             sPL.domain    = @(coor) isTop(coor);
%             sPL.direction = 2;
%             sPL.value     = 200;
%             s.pointloadFun = DistributedLoad(obj.mesh, sPL);%DistributedLoad(obj.mesh, sPL);

%             [bM,l2g] = obj.mesh.getBoundarySubmesh(sPL.domain);
% 
%             sAF.fHandle = @(x) [0*x(1,:,:);sPL.value*ones(size(x(1,:,:)))];
%             sAF.ndimf   = 2;
%             sAF.mesh    = bM;
%             xFun = AnalyticalFunction(sAF);
%             xFunP1  =xFun.project('P1');
% 
%             s.mesh = bM;
%             s.type = 'ShapeFunction';
%             s.quadType = 2;
%             rhsI       = RHSintegrator.create(s);
%             test = LagrangianFunction.create(bM,xFun.ndimf,'P1');
%             Fext2 = rhsI.compute(xFunP1,test);   
%             Fext3 = reshape(Fext2,[bM.ndim,bM.nnodes])';
% 
%             Fext = zeros(obj.mesh.nnodes,2);
%             Fext(l2g,:) = Fext3;
% 
% %             obj.FextInitial = Fext; 
%             [pl_dofs, ~, pl_vals] = find(Fext);
%             ndofs = obj.mesh.nnodes*obj.mesh.ndim;
%             fVals = zeros(ndofs,1);
%             fVals(pl_dofs) = pl_vals;
%             fVals = obj.reshapeToMatrix(fVals);

%             s.pointloadFun.fun = LagrangianFunction.create(obj.mesh,xFun.ndimf,'P1');
%             s.pointloadFun.fun.fValues = fVals;
%             s.pointloadFun.fValues = fVals;
%             s.pointloadFun.values = pl_vals;
            s.periodicFun  = [];
            s.mesh         = obj.mesh;

            bc = BoundaryConditions(s);
%             obj.boundaryConditions = bc;
        end

        function rshp = reshapeToMatrix(obj, A)
            rshp = reshape(A,[obj.mesh.ndim,obj.mesh.nnodes])';
        end

    end

end