classdef Tutorial02p3ElasticityAMG < handle

% WARNING!
% - Install a valid python version (3.8 for example)
% - Install pyAMG python library
% - Install further python dependencies (numpy, scipy, ...)

    properties (Access = private)
        mesh
        young
        poisson
        material
        stateProblem
    end

    methods (Access = public)

        function obj = Tutorial02p3ElasticityAMG()
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
            obj.mesh = UnitTriangleMesh(50,50);
        end

        function computeElasticProperties(obj)
            E  = 1;
            nu = 1/3;
            obj.young   = ConstantFunction.create(E,obj.mesh);
            obj.poisson = ConstantFunction.create(nu,obj.mesh);
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
            s.solverCase = obj.createSolver(s);
            fem = ElasticProblem(s);
            fem.solve();
            obj.stateProblem = fem;
        end

        function solver = createSolver(obj,s)
            BCAp = BCApplier(s);
            Rfull  = obj.computeRigidBodyModes([0.5,0.5]);
            for i = 1:size(Rfull,2)
                R(:,i) = BCAp.fullToReducedVectorDirichlet(Rfull(:,i));
            end
            s.type = 'ELASTIC';
            s.nullSpace = R;
            s.nLevels = 5;
            s.tol = 1e-8;
            s.maxIter = 1;
            p     = pyAMG.create(s);

            sS.preconditioner = p;
            sS.tol = 1e-5;
            solver = PCG(sS);
        end

        function R = computeRigidBodyModes(obj,refPoint)
            rigModes = RigidBodyFunction.create(obj.mesh,refPoint);
            RFun = rigModes.projectBasisFunctions('P1');
            for i = 1:length(RFun)
                R(:,i) = reshape(RFun{i}.fValues',[],1);
            end
        end

        function bc = createBoundaryConditions(obj)
            xMax    = max(obj.mesh.coord(:,1));
            yMax    = max(obj.mesh.coord(:,2));
            isDir   = @(coor)  abs(coor(:,1))==0;
            isForce = @(coor)  (abs(coor(:,1))==xMax & abs(coor(:,2))>=0.3*yMax & abs(coor(:,2))<=0.7*yMax);

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
                pl = TractionLoad(obj.mesh, sPL{i}, 'DIRAC');
                pointloadFun = [pointloadFun, pl];
            end
            s.pointloadFun = pointloadFun;

            s.periodicFun  = [];
            s.mesh = obj.mesh;
            bc = BoundaryConditions(s);
        end

    end

end