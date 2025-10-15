classdef TopOptTests < handle & matlab.unittest.TestCase

% En el futur: agrupar testsTO malla en comú i sense definir functionals, opt..;
% crear difs funcions a topopttests depenent tipus functionals...

    properties (TestParameter)
        testsTO = {
            'test_cantilever', 'test_cantilever2', 'test_cantilever3', ...
            'testDualNestedInPrimal_WithProjectedGradient', ...
            'testDualNestedInPrimal_WithSlerp', ...
            'test_interiorPerimeter',...
            'test_anisotropy','test_anisotropy_interior','test_nullspace',...
            'test_interiorPerimeterPDErho','test_filterLump','test_cantilever_IPM',...
            'test_dirichletProjection','test_gripping','test_micro', 'test_micro2',...
            'test_micro3','test_boundFormFilterAndProject','test_cantilever_SIMPP3',...
            'test_infillLS','test_isoPerLS', 'test_PerpnormLS','test_segment','test_droplet'
            }
    end

    methods (Test, TestTags = {'TopOpt', 'Various', 'testsTO'})

        function testFastDisplacement(testCase, testsTO)
            t = testCase.runTopOptTest(testCase,testsTO);
            xNew = t.designVariable.fun.fValues;
            load([testsTO,'.mat'],'x');
            err = norm(x - xNew)/norm(x);
            tol = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol);
        end

    end

    methods (Static, Access = private)
        function t = runTopOptTest(obj,testName)
            run(testName);
            m      = gid.mesh;
            dim    = gid.dim;
            bc     = gid.boundaryConditions;
            x      = obj.createDesignVariable(designVariable,m,geomFunSettings,plotting);
            filtersCost = obj.createFilters(filterCostType,m,filterCostSettings);
            filtersConstraint = obj.createFilters(filterConstraintType,m,filterConstraintSettings);
            mI     = obj.createMaterialInterpolator(materialType,method,m,dim);
            mat    = obj.createMaterial(x,mI,m);
            fem    = obj.createElasticProblem(m,mat,ptype,dim,bc);
            Msmooth = obj.createMassMatrix(m,x);
            if exist('micro')
                s = micro;
            else
                s = [];
            end
            base = obj.createChiOmega(m);
            sFCost = obj.createCost(cost,weights,m,fem,filtersCost,mat,Msmooth,filename,s,base);
            sFConstraint = obj.createConstraint(constraint,target,m,fem,filtersConstraint,mat,Msmooth,base);
            l.nConstraints = length(constraint);
            lam    = DualVariable(l);
            primal = optimizerUnconstrained;
            primalUp = obj.createPrimalUpdater(m,x,primal);
            t      = obj.createOptimizer(optimizer,primalUp,primal,monitoring,sFCost,sFConstraint,x,lam,maxiter,constraint_case,target);
            close all;
        end

        function x = createDesignVariable(type,mesh,gSet,plotting)
            g      = GeometricalFunction(gSet);
            lsFun  = g.computeLevelSetFunction(mesh);
            switch type
                case {'Density','DensityAndBound'}
                    ss.fValues = 1 - heaviside(lsFun.fValues);
                    ss.mesh    = mesh;
                    ss.order   = 'P1';
                    fun        = LagrangianFunction(ss);
                case 'LevelSet'
                    fun = lsFun;
            end
            s.fun  = fun;
            s.mesh = mesh;
            s.type = type;
            s.plotting = plotting;
            x      = DesignVariable.create(s);
        end

        function filtersCost = createFilters(type,mesh,settings)
            for i = 1:length(type)
                s            = settings{i};
                s.filterType = type{i};
                s.mesh       = mesh;
                s.test  = LagrangianFunction.create(mesh,1,'P0');
                s.trial = LagrangianFunction.create(mesh,1,'P1');
                if isempty(type{i})
                    filtersCost{i} = [];
                else
                    filtersCost{i} = Filter.create(s);
                end
            end
        end

        function mI = createMaterialInterpolator(materialType,method,mesh,dim)
            E1   = 1;
            E0   = 1e-3;
            nu1  = 1/3;
            nu0  = 1/3;
            ndim = mesh.ndim;

            matA.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E0,nu0);
            matA.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E0,nu0,ndim);

            matB.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E1,nu1);
            matB.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E1,nu1,ndim);

            s.typeOfMaterial = materialType;
            s.interpolation  = method;
            s.dim            = dim;
            s.matA = matA;
            s.matB = matB;

            mI = MaterialInterpolator.create(s);
        end

        function m = createMaterial(x,mI,m)
            switch class(x)
                case 'DensityAndBound'
                    f = x.density.obtainDomainFunction();
                otherwise
                    f = x.obtainDomainFunction();
            end
            f = f{1}.project('P1');
            s.type                 = 'DensityBased';
            s.density              = f;
            s.materialInterpolator = mI;
            s.dim                  = '2D';
            s.mesh                 = m;
            m = Material.create(s);
        end

        function fem = createElasticProblem(mesh,mat,scale,dim,bc)
            s.mesh               = mesh;
            s.scale              = scale;
            s.material           = mat;
            s.dim                = dim;
            s.boundaryConditions = bc;
            s.interpolationType  = 'LINEAR';
            s.solverType         = 'REDUCED';
            switch s.scale
                case 'MACRO'
                    s.solverMode = 'DISP';
                case 'MICRO'
                    s.solverMode = 'FLUC';
            end
            s.solverCase         = DirectSolver();
            s.type               = 'ELASTIC';
            fem                  = PhysicalProblem.create(s);
        end

        function M = createMassMatrix(mesh,x)
            vF = LagrangianFunction.create(mesh,1,'P1');
            M = IntegrateLHS(@(u,v) DP(v,u),vF,vF,mesh,'Domain'); 
            M = eye(size(M));
            switch class(x)
                case 'DensityAndBound'
                    n = mesh.nnodes;
                    M(n+1,n+1) = 1;
            end
        end

        function base = createChiOmega(mesh)
            levelSet         = -ones(mesh.nnodes,1);
            s.backgroundMesh = mesh;
            s.boundaryMesh   = mesh.createBoundaryMesh();
            base             = UnfittedMesh(s);
            base.compute(levelSet);
        end

        function sFCost = createCost(cost,weights,mesh,fem,filter,mat,Msmooth,filename,s,base)
            for i = 1:length(cost)
                s.type            = cost{i};
                s.mesh            = mesh;
                s.physicalProblem = fem;
                s.filter          = filter{i};
                s.material        = mat;
                s.filename        = filename;
                s.base            = base;
                sF{i}             = ShapeFunctional.create(s);
            end
            ss.shapeFunctions = sF;
            ss.weights        = weights;
            ss.Msmooth        = Msmooth;
            sFCost            = Cost(ss);
        end

        function sFConstraint = createConstraint(constraint,target,mesh,fem,filter,mat,Msmooth,base)
            k = 1;
            for i = 1:length(constraint)
                switch class(filter{k})
                    case 'FilterAndProject'
                        s.filterDesignVariable = filter{k};
                        s.filterGradient       = filter{k+1};
                        k = k+2;
                    otherwise
                        s.filter = filter{k};
                        k = k+1;
                end
                s.type            = constraint{i};
                s.target          = target(i);
                s.mesh            = mesh;
                s.physicalProblem = fem;
                s.material        = mat;
                s.gradientTest    = LagrangianFunction.create(mesh,1,'P1');
                s.base            = base;
                sF{i}             = ShapeFunctional.create(s);
            end
            ss.shapeFunctions = sF;
            ss.Msmooth        = Msmooth;
            sFConstraint      = Constraint(ss);
        end

        function p = createPrimalUpdater(m,x,primal)
            if isempty(primal)
                p = [];
            else
                s.mesh           = m;
                s.designVariable = x;
                s.primal         = primal;
                switch x.type
                    case 'Density'
                        s.ub = 1;
                        s.lb = 0;
                    case 'DensityAndBound'
                        s.ub = [ones(x.density.fun.mesh.nnodes,1);1000];
                        s.lb = [zeros(x.density.fun.mesh.nnodes,1);-1000];
                end
                s.tauMax = 1000;
                s.tau    = [];
                p        = PrimalUpdaterFactory.create(s);
            end
        end

        function s = createOptimizer(type,primalUp,primal,monitoring,cost,constraint,x,lam,maxIter,constraintCase,target)
            s.type           = type;
            s.monitoring     = monitoring;
            s.cost           = cost;
            s.constraint     = constraint;
            s.designVariable = x;
            s.dualVariable   = lam;
            s.maxIter        = maxIter;
            s.tolerance      = 1e-8;
            s.constraintCase = constraintCase;
            s.volumeTarget   = target; % will dissappear
            s.primalUpdater  = primalUp;
            s.primal         = primal;
            s.etaNorm        = 0.05;
            s.gJFlowRatio    = 1; % Only NullSpace
            s.tauMax         = 1000;
            s.etaMax         = 1;
            s.etaMaxMin      = 0.01;
            s.etaNormMin     = 0.05;
            switch x.type
                case 'Density'
                    s.ub = 1;
                    s.lb = 0;
                case 'DensityAndBound'
                    s.ub = [ones(x.density.fun.mesh.nnodes,1);1000];
                    s.lb = [zeros(x.density.fun.mesh.nnodes,1);-1000];
            end
            opt = Optimizer.create(s);
            opt.solveProblem();
        end
    end
end