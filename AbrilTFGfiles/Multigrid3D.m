classdef Multigrid3D < handle

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
        nSubdomains
        tolSameNode
    end

    methods (Access = public)

        function obj = Multigrid3D()
            obj.init();
            obj.createMesh();
            obj.computeElasticProperties();
            obj.createMaterial();
            obj.solveElasticProblem();
        end

    end

    methods (Access = private)

        function init(obj)
            obj.nSubdomains = [35 1 1];
            obj.tolSameNode = 1e-9;
        end

        function createMesh(obj)    
            obj.mesh = UnitTriangleMesh(50,50);
            filename = 'DEF_Q8_wing_1.mat';
            load(filename);
            s.coord    = EIFEoper.MESH.COOR;
            s.connec   = EIFEoper.MESH.CN;

            maxC= max(s.coord);
            minC = min(s.coord);

            delta=1E-5;  % 1E-9--> general case   1E-5 --> airfoil

            s.coord(s.coord(:,1)== maxC(1) & s.coord(:,3)==maxC(3),:) =...
                s.coord(s.coord(:,1)== maxC(1) & s.coord(:,3)==maxC(3),:)-[0,0,delta];

            s.coord(s.coord(:,1)== maxC(1) & s.coord(:,3)==minC(3),:) =...
                s.coord(s.coord(:,1)== maxC(1) & s.coord(:,3)==minC(3),:)+[0,0,delta];

            s.coord(s.coord(:,1)== minC(1) & s.coord(:,3)==maxC(3),:) =...
                s.coord(s.coord(:,1)== minC(1) & s.coord(:,3)==maxC(3),:)-[0,0,delta];

            s.coord(s.coord(:,1)== minC(1) & s.coord(:,3)==minC(3),:) =...
                s.coord(s.coord(:,1)== minC(1) & s.coord(:,3)==minC(3),:)+[0,0,delta];

            s.coord(s.coord(:,2)== maxC(2) & s.coord(:,3)==maxC(3),:) =...
                s.coord(s.coord(:,2)== maxC(2) & s.coord(:,3)==maxC(3),:)-[0,0,delta];

            s.coord(s.coord(:,2)== maxC(2) & s.coord(:,3)==minC(3),:) =...
                s.coord(s.coord(:,2)== maxC(2) & s.coord(:,3)==minC(3),:)+[0,0,delta];

            s.coord(s.coord(:,2)== minC(2) & s.coord(:,3)==maxC(3),:) =...
                s.coord(s.coord(:,2)== minC(2) & s.coord(:,3)==maxC(3),:)-[0,0,delta];

            s.coord(s.coord(:,2)== minC(2) & s.coord(:,3)==minC(3),:) =...
                s.coord(s.coord(:,2)== minC(2) & s.coord(:,3)==minC(3),:)+[0,0,delta];

            mS       = Mesh.create(s);
            mD       = obj.createMeshDomain(mS);
            obj.mesh = mD;

        end

        function mD= createMeshDomain(obj,mR)
            s.nsubdomains   = obj.nSubdomains; %  num along x and y
            s.meshReference = mR;
            s.tolSameNode = obj.tolSameNode;
            m = MeshCreatorFromRVE3D(s);
            [mD,mSb,iC,~,lG,iCR,discmesh] = m.create();
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
            s.boundaryConditions = obj.createBoundaryConditions(obj.mesh);
            s.solverType = 'REDUCED';
            s.solverMode = 'DISP';
            s.solverCase = obj.createSolver(s);
            fem = ElasticProblem(s);
            fem.solve();
            obj.stateProblem = fem;
        end

        function solver = createSolver(obj,s)
            BCAp = BCApplier(s);
            % Rfull  = obj.computeRigidBodyModes([-0.05,0.0011,0.0152]);
            Rfull  = obj.computeRigidBodyModes([-0.05,0.0011,17]);
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
            sS.tol = 1e-8;
            solver = PCG2(sS);
        end

        function R = computeRigidBodyModes(obj,refPoint)
            rigModes = RigidBodyFunction.create(obj.mesh,refPoint);
            RFun = rigModes.projectBasisFunctions('P1');
            for i = 1:length(RFun)
                R(:,i) = reshape(RFun{i}.fValues',[],1);
            end
        end


         function [bc,Dir,PL] = createBoundaryConditions(obj,mesh)
            [Dir,PL]  = obj.createRawBoundaryConditions();
            dirichletFun = [];
            for i = 1:numel(Dir)
                dir = DirichletCondition(obj.mesh, Dir{i});
                dirichletFun = [dirichletFun, dir];
            end
            pointload = TractionLoad(mesh,PL,'DIRAC');
            s.pointloadFun = pointload;
            s.dirichletFun = dirichletFun;
            s.periodicFun  =[];
            s.mesh         = mesh;
            bc             = BoundaryConditions(s);
         end

         function [Dir,PL] = createRawBoundaryConditions(obj)
            minx = min(obj.mesh.coord(:,1));
            maxx = max(obj.mesh.coord(:,1));
            miny = min(obj.mesh.coord(:,2));
            maxy = max(obj.mesh.coord(:,2));
            minz = min(obj.mesh.coord(:,3));
            maxz = max(obj.mesh.coord(:,3));
            tolBound = obj.tolSameNode;
            isLeft   = @(coor) (abs(coor(:,1) - minx)   < tolBound);
            isRight  = @(coor) (abs(coor(:,1) - maxx)   < tolBound);
            isBottom = @(coor) (abs(coor(:,3) - minz)   < tolBound);
            
            Dir{1}.domain    = @(coor) isLeft(coor);%| isRight(coor) ;
            Dir{1}.direction = [1,2,3];
            Dir{1}.value     = 0;

            PL.domain    = @(coor) isRight(coor);
            PL.direction = 2;        % 3--> general    2--> Airfoil
            PL.value     = -1;       %Set displacement intensity 
        end 


        

    end

end