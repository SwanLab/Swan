classdef MultimaterialTesting < handle

    properties (Access = private)
        mesh
        designVariable
        filter
        physicalProblem
        compliance
        volumeA
        volumeB
        volumeC
        cost
        constraint
        dualVariable
        optimizer
        nMat
        pdeCoeff
        matInterp
        mat
        bc
        nLevelSet
    end

    methods (Access = public)

        function obj = MultimaterialTesting()
            obj.init()
            obj.createMesh();
            obj.createDesignVariable();
            obj.createMaterialProperties();
            obj.createInterpolators();
            obj.createBoundaryConditions();
            obj.createElasticProblem();
            obj.createCompliance();
            obj.createVolumeConstraints();
            obj.createCost();
            obj.createConstraint();
            obj.createDualVariable();
            obj.createOptimizer();
        end

        function lsVals = solve(obj)
            obj.optimizer.solveProblem();
            lsVals = obj.designVariable.fun.fValues;
        end

    end

    methods (Access = private)

        function init(obj)
            close all;
            obj.nMat = 4; % including weakest mat
            obj.nLevelSet = 3;
        end

        function createMesh(obj)
            obj.mesh = TriangleMesh(2,1,20,10);
        end

        function createDesignVariable(obj)
            s.mesh                 = obj.mesh;
            s.type                 = 'Full';
            lsFun{1}               = @(x) -ones(size(x(1,:,:)));
            lsFun{2}               = @(x) -cos(x(1,:,:))+0.5;
            lsFun{3}               = @(x) sin(x(1,:,:))-0.5;
            
            s.type                 = 'MultiLevelSet';
            s.lsFun                = lsFun;
            s.mesh                 = obj.mesh;
            s.unitM                = obj.createMassMatrix();
            s.plotting             = true;
            obj.designVariable     = DesignVariable.create(s);
        end

        function matI = createSingleMaterialProperty(obj,E,nu)
            mu          = E./(2.*(1+nu));
            la          = nu.*E./((1+nu).*(1-2.*nu)); % plain strain
            matI.young  = E;
            matI.nu     = nu;
            matI.shear  = mu;
            matI.lambda = 2.*mu.*la./(la+2.*mu); % plane stress
        end

        function createMaterialProperties(obj)
            obj.mat.A = obj.createSingleMaterialProperty(200E9,0.25);
            obj.mat.B = obj.createSingleMaterialProperty(100E9,0.25);
            obj.mat.C = obj.createSingleMaterialProperty(50E9,0.25);
            obj.mat.D = obj.createSingleMaterialProperty(0.2E9,0.25);
        end
        
        function createInterpolators(obj)
            s.mat        = obj.mat;
            s.m          = obj.mesh;
            
            obj.pdeCoeff = PDECoefficientsComputer(s);  


            sC.E  = [obj.mat.A.young,obj.mat.B.young,obj.mat.C.young,obj.mat.D.young];
            sC.CA = obj.pdeCoeff.tensor(:,1);

            E   = AnalyticalFunction.create(@(x) 100E9*ones(size(squeezeParticular(x(1,:,:),2))),1,obj.mesh);
            nu  = AnalyticalFunction.create(@(x) 0.25*ones(size(squeezeParticular(x(1,:,:),2))),1,obj.mesh);
            s.type    = 'ISOTROPIC';
            s.ptype   = 'ELASTIC';
            s.ndim    = obj.mesh.ndim;
            s.young   = E;
            s.poisson = nu;
            tensor    = Material.create(s);
            tensorEv  = tensor.evaluate([0;0]);
            CB     = tensorEv(:,:,1,1);
            sC.CB  = zeros(size(CB));
            sC.CB(2,2) = 2*CB(3,3);
            sC.CB(1,1) = CB(1,1);
            sC.CB(3,3) = CB(2,2);
            sC.CB(1,3) = CB(1,2);
            sC.CB(3,1) = CB(2,1);
            obj.matInterp = MultiMaterialInterpolation(sC);
        end

        function createBoundaryConditions(obj)
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
            obj.bc = BoundaryConditions(s);
        end

        function createElasticProblem(obj)
            s.mesh = obj.mesh;
            s.scale = 'MACRO';
            s.material = obj.createMaterial();
            s.dim = '2D';
            s.boundaryConditions = obj.bc;
            %s.interpolationType = 'LINEAR';
            s.solverType = 'REDUCED';
            s.solverMode = 'DISP';
            s.solverCase = 'DIRECT';
            fem = ElasticProblem(s);
            obj.physicalProblem = fem;
        end

        function createCompliance(obj)
            s.nMat = obj.nMat;
            s.mesh = obj.mesh;
            s.stateProblem = obj.physicalProblem;
            s.material = obj.createMaterial();
            s.materialInterpolator = obj.matInterp;
            c = MultiMaterialComplianceFunctional(s);
            obj.compliance = c;
        end

        function createVolumeConstraints(obj)
            obj.volumeA = obj.createIndivVolumeConstraint(0.1,1);
            obj.volumeB = obj.createIndivVolumeConstraint(0.1,2);
            obj.volumeC = obj.createIndivVolumeConstraint(0.1,3);
        end

        function v = createIndivVolumeConstraint(obj,target,ID)
            s.volumeTarget = target;
            s.nMat         = 4;
            s.matID        = ID;
            s.mesh         = obj.mesh;
            v              = MultiMaterialVolumeConstraint(s);
         end

         function createCost(obj)
            s.shapeFunctions{1} = obj.compliance;
            s.weights           = 1;
            s.Msmooth           = obj.createMassMatrix();
            obj.cost            = Cost(s);
         end

         function createConstraint(obj)
            s.shapeFunctions{1} = obj.volumeA;
            s.shapeFunctions{2} = obj.volumeB;
            s.shapeFunctions{3} = obj.volumeC;
            s.Msmooth           = obj.createMassMatrix();
            obj.constraint      = Constraint(s);
         end

         function createDualVariable(obj)
            s.nConstraints   = obj.nMat-1;
            l                = DualVariable(s);
            obj.dualVariable = l;
         end

         function createOptimizer(obj)
            s.monitoring     = true;
            s.cost           = obj.cost;
            s.constraint     = obj.constraint;
            s.designVariable = obj.designVariable;
            s.dualVariable   = obj.dualVariable;
            s.maxIter        = 2;
            s.tolerance      = 1e-8;
            s.constraintCase = repmat({'EQUALITY'},[obj.nMat,1]);
            s.primal         = 'SLERP';
            s.ub             = inf;
            s.lb             = -inf;
            s.etaNorm        = inf;
            s.gJFlowRatio    = 2;
            obj.optimizer    = OptimizerNullSpace(s);
        end

        function m = createMaterial(obj)
            x                      = obj.designVariable;
            s.type                 = 'DensityBased';
            s.density              = x;
            s.materialInterpolator = obj.matInterp;
            s.dim                  = '2D';
            m = Material.create(s);
        end

        function M = createMassMatrix(obj)
            nnodes  = obj.mesh.nnodes*obj.nLevelSet;
            %nnodes  = obj.mesh.nnodes;
            indices = transpose(1:nnodes);
            vals    = ones(size(indices));
            h       = obj.mesh.computeMeanCellSize();
            M       = h^2*sparse(indices,indices,vals,nnodes,nnodes);
        end


    end
end