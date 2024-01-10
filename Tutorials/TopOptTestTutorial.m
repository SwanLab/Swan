classdef TopOptTestTutorial < handle

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

        function obj = TopOptTestTutorial()
            obj.init()
            obj.createMesh();
            obj.createDesignVariable();            
            obj.createFilter();
            obj.createMaterialInterpolator();
            obj.createElasticProblem();
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
            s.fHandle = @(x) ones(size(squeezeParticular(x(1,:,:),1)));
            s.ndimf   = 1;
            s.mesh    = obj.mesh;
            aFun      = AnalyticalFunction(s);            
            s.fun     = aFun.project('P1');
            s.mesh    = obj.mesh;                        
            s.type = 'Density';
            dens    = DesignVariable.create(s);   
            obj.designVariable = dens;
        end

        function createFilter(obj)
            s.filterType = 'LUMP';
            s.mesh  = obj.mesh;
            s.trial = P1Function.create(obj.mesh,1);
            f = Filter.create(s);
            obj.filter = f;
        end       

        function createMaterialInterpolator(obj)
            E0 = 1e-3;
            nu0 = 1/3;
            E = AnalyticalFunction.create(@(x) E0*ones(size(squeeze(x(1,:,:)))),1,obj.mesh);
            nu = AnalyticalFunction.create(@(x)nu0*ones(size(squeeze(x(1,:,:)))),1,obj.mesh);
            tensorA = obj.createMaterial(E,nu);            
             
            E1  = 1;
            nu1 = 1/3;            
            E  = AnalyticalFunction.create(@(x)E1*ones(size(squeeze(x(1,:,:)))),1,obj.mesh);
            nu = AnalyticalFunction.create(@(x)nu1*ones(size(squeeze(x(1,:,:)))),1,obj.mesh);
            tensorB = obj.createMaterial(E,nu);

            ndim = 2;
            matA.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E0,nu0);
            matA.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E0,nu0,ndim);

            matB.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E1,nu1);
            matB.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E1,nu1,ndim);

            s.typeOfMaterial = 'ISOTROPIC';
            s.interpolation  = 'SIMPALL';
            s.dim            = '2D';
            s.matA = matA;
            s.matB = matB;

            m = MaterialInterpolator.create(s);
            obj.materialInterpolator = m;            
        end    

        function m = createMaterial(obj,young,poisson)
            s.type    = 'ISOTROPIC';
            s.ptype   = 'ELASTIC';
            s.ndim    = obj.mesh.ndim;
            s.young   = young;
            s.poisson = poisson;
            m = Material.create(s);               
        end


        function createElasticProblem(obj)
            s.mesh = obj.mesh;
            s.scale = 'MACRO';
            s.material = obj.createInterpolatedMaterial(obj.designVariable.fun);
            s.dim = '2D';
            s.bc = obj.createBoundaryConditions();
            s.interpolationType = 'LINEAR';
            fem = ElasticProblem(s);
            obj.physicalProblem = fem;
        end

        function createCompliance(obj)
            s.mesh                 = obj.mesh;
            s.filter               = obj.filter;
            s.stateProblem         = obj.physicalProblem;
            s.materialInterpolator = obj.materialInterpolator;
            c                      = ComplianceFunctional(s);
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
            s.maxIter        = 100;
            s.tolerance      = 1e-8;
            s.constraintCase = 'EQUALITY';
            s.ub             = 1;
            s.lb             = 0;
            opt = OptimizerMMA(s);
            opt.solveProblem();
            obj.optimizer = opt;
        end

        function mat = createInterpolatedMaterial(obj,dens)
            mI   = obj.materialInterpolator;
            mat  = mI.compute(dens);
        end
        
        function bc = createBoundaryConditions(obj)
            bM = obj.mesh.createBoundaryMesh();

            dBC.boundaryId   = 1;
            dBC.dof          = [1,2];
            dBC.value        = [0,0];
            nBC.boundaryId   = 2;
            nBC.dof          = 2;
            nBC.value        = -1;

            [dirichlet,pointload] = obj.createBc(bM,dBC,nBC);
            bc.dirichlet=dirichlet;
            bc.pointload=pointload;
        end

       function [dirichlet,pointload] = createBc(obj,bMesh,dBC,nBC)
            dirichlet = obj.createBondaryCondition(bMesh,dBC);
            pointload = obj.createBondaryCondition(bMesh,nBC);
        end

        function cond = createBondaryCondition(obj,bM,condition)
            nbound = length(condition.boundaryId);
            cond = zeros(1,3);
            for ibound=1:nbound
                ncond  = length(condition.dof(nbound,:));
                nodeId = unique(bM{condition.boundaryId(ibound)}.globalConnec);
                nbd   = length(nodeId);
                for icond=1:ncond
                    bdcond= [nodeId, repmat(condition.dof(icond),[nbd,1]), repmat(condition.value(icond),[nbd,1])];
                    cond=[cond;bdcond];
                end
            end
            cond = cond(2:end,:);
        end        

    end

end