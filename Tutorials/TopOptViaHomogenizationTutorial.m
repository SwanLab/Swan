classdef TopOptViaHomogenizationTutorial < handle

    properties (Access = private)
        mesh
        filter
        designVariable
        materialInterpolator
        physicalProblem
        compliance
        volume
        cost
        constraint
        dualVariable
        optimizer
    end

    methods (Access = public)

        function obj = TopOptViaHomogenizationTutorial()
            obj.init()
            obj.createMesh();
            obj.createDesignVariable();            
            obj.createFilter();
            obj.createMaterialInterpolator();
            obj.createElasticProblem();
            obj.createComplianceFromConstiutive();            
            obj.createCompliance();
            obj.createVolume();
            obj.createCost();
            obj.createConstraint();
            obj.createDualVariable();
            obj.createOptimizer();
        end

    end

    methods (Access = private)

        function init(obj)

        end

        function createMesh(obj)
            %UnitMesh better
            x1      = linspace(0,1,50);
            x2      = linspace(0,1,50);
            [xv,yv] = meshgrid(x1,x2);
            [F,V]   = mesh2tri(xv,yv,zeros(size(xv)),'x');
            s.coord  = V(:,1:2);
            s.connec = F;
            obj.mesh = Mesh(s);            
        end

        function createDesignVariable(obj)
            s.fHandle = @(x) 0.5*ones(size(squeezeParticular(x(1,:,:),1)));
            s.ndimf   = 1;
            s.mesh    = obj.mesh;
            aFun      = AnalyticalFunction(s);            
            s.fun{1}  = aFun.project('P1');
            s.fun{2}  = aFun.project('P1');
            s.mesh    = obj.mesh;                        
            s.type    = 'MicroParams';
            desVar    = DesignVariable.create(s);   
            obj.designVariable = desVar;
        end

        function createFilter(obj)
            s.filterType = 'LUMP';
            s.mesh  = obj.mesh;
            s.trial = P1Function.create(obj.mesh,1);
            f = Filter.create(s);
            obj.filter = f;
        end       

        function createMaterialInterpolator(obj)
%             ndim = 2;            
%             E0 = 1e-3; 
%             nu0 = 1/3;
%             matA.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E0,nu0);
%             matA.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E0,nu0,ndim);
% 
%             E1 = 1;
%             nu1 = 1/3;              
%             matB.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E1,nu1);
%             matB.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E1,nu1,ndim);

            s.interpolation  = 'HomogenizedMicrostructure';
            s.fileName = 'Rectangle';

            m = MaterialInterpolator.create(s);
            obj.materialInterpolator = m;            
        end    

        function createElasticProblem(obj)
            s.mesh = obj.mesh;
            s.scale = 'MACRO';
            s.material = obj.createInterpolatedMaterial(obj.designVariable.fun);
            s.dim = '2D';
            s.bc = obj.createBoundaryConditions();
            s.solverType = 'REDUCED';
            s.solverMode = 'DISP';
            s.interpolationType = 'LINEAR';
            fem = ElasticProblem(s);
            obj.physicalProblem = fem;
        end

        function createComplianceFromConstiutive(obj)
            s.mesh         = obj.mesh;
            s.stateProblem = obj.physicalProblem;
            c = ComplianceFromConstiutiveTensor(s);
            obj.compliance = c;
        end        

        function createCompliance(obj)
            s.mesh                 = obj.mesh;
            s.filter               = obj.filter;
            s.stateProblem         = obj.physicalProblem;
            s.materialInterpolator = obj.materialInterpolator;
            c                      = ComplianceFunctionalFromVademecum(s);
            obj.compliance = c;
        end

        function createVolume(obj)
            s.mesh   = obj.mesh;
            s.filter = obj.filter;
            s.volumeTarget = 0.4;
            v = VolumeConstraint(s);
            obj.volume = v;
        end

        function createCost(obj)
            s.ndof              = obj.mesh.nnodes;
            s.shapeFunctions{1} = obj.compliance;
            s.weights           = 1;
            obj.cost            = Cost(s);
        end

        function createConstraint(obj)
            s.ndof              = obj.mesh.nnodes;
            s.shapeFunctions{1} = obj.volume;
            obj.constraint      = Constraint(s);
        end

        function createDualVariable(obj)
            s.nConstraints   = 1;
            l                = DualVariable(s);
            obj.dualVariable = l;
        end

        function createOptimizer(obj)
            s.cost           = obj.cost;
            s.constraint     = obj.constraint;
            s.designVariable = obj.designVariable;
            s.dualVariable   = obj.dualVariable;
            s.maxIter        = 1000;
            s.tolerance      = 1e-8;
            s.constraintCase = 'EQUALITY';
            s.ub             = 1;
            s.lb             = 0;
            opt = OptimizerMMA(s);
            opt.solveProblem();
            obj.optimizer = opt;
        end

        function mat = createInterpolatedMaterial(obj,desVar)
            mI   = obj.materialInterpolator;
            mat  = mI.computeConsitutiveTensor(desVar);
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
            bc.dirichletFun = dirichletFun;

            pointloadFun = [];
            for i = 1:numel(sPL)
                pl = PointLoad(obj.mesh, sPL{i});
                pointloadFun = [pointloadFun, pl];
            end
            bc.pointloadFun = pointloadFun;

            bc.periodicFun  = [];
        end    

    end

end