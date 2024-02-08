classdef TopOptTests < handle & matlab.unittest.TestCase

    properties (TestParameter)
        testsTO = {
            'test_cantilever', 'test_cantilever2', 'test_cantilever3', ...
            'test_micro', 'test_micro2', ...
            'testDualNestedInPrimal_WithProjectedGradient', ...
            'testDualNestedInPrimal_WithSlerp', ...
            'test_interiorPerimeter',...
            'test_anisotropy','test_anisotropy_interior','test_nullspace','test_gripping',...
            'test_interiorPerimeterPDErho','test_filterLump','test_cantilever_IPM'
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
            gid    = obj.readGidFile(filename);
            m      = gid.mesh;
            dim    = gid.dim;
            bc     = gid.bc;
            x      = obj.createDesignVariable(initialCase,designVariable,m);
            filtersCost = obj.createFilters(filterCostType,m);
            filtersConstraint = obj.createFilters(filterConstraintType,m);
            mI     = obj.createMaterialInterpolator(materialType,method,m,E1,E0,nu1,nu0,dim);
            fem    = obj.createElasticProblem(x,m,mI,ptype,dim,bc);
            sFCost = obj.createCost(cost,weights,m,fem,filtersCost,mI);
            sFConstraint = obj.createConstraint(constraint,target,m,fem,filtersConstraint,mI);
            l.nConstraints = length(constraint);
            lam    = DualVariable(l);
            primal = optimizerUnconstrained;
            t      = obj.createOptimizer(optimizer,primal,monitoring,sFCost,sFConstraint,x,lam,maxiter,constraint_case,target);
        end

        function s = readGidFile(file)
            a.fileName = file;
            s          = FemDataContainer(a);
        end

        function x = createDesignVariable(initCase,type,mesh)
            s.type = initCase;
            g      = GeometricalFunction(s);
            lsFun  = g.computeLevelSetFunction(mesh);
            s.fun  = lsFun;
            s.mesh = mesh;
            s.type = type;
            x      = DesignVariable.create(s);
        end

        function filtersCost = createFilters(type,mesh)
            for i = 1:length(type)
                s.filterType = type{i};
                s.mesh       = mesh;
                s.test       = P0Function.create(mesh,1);
                s.trial      = P1Function.create(mesh,1);
                if isempty(type{i})
                    filtersCost{i} = [];
                else
                    filtersCost{i} = Filter.create(s);
                end
            end
        end

        function mI = createMaterialInterpolator(materialType,method,mesh,E1,E0,nu1,nu0,dim)
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

        function fem = createElasticProblem(x,mesh,mI,scale,dim,bc)
            f   = x.obtainDomainFunction();
            f   = f.project('P1');
            mat = mI.computeConsitutiveTensor(f);
            s.mesh              = mesh;
            s.scale             = scale;
            s.material          = mat;
            s.dim               = dim;
            s.bc                = bc;
            s.interpolationType = 'LINEAR';
            s.solverType        = 'REDUCED';
            s.solverMode        = 'DISP';
            fem                 = ElasticProblem(s);
        end

        function sFCost = createCost(cost,weights,mesh,fem,filter,mI)
            for i = 1:length(cost)
                s.type                 = cost{i};
                s.mesh                 = mesh;
                s.physicalProblem      = fem;
                s.filter               = filter{i};
                s.materialInterpolator = mI;
                sF{i}                  = ShapeFunctional.create(s);
            end
            ss.shapeFunctions = sF;
            ss.weights        = weights;
            sFCost            = Cost(ss);
        end

        function sFConstraint = createConstraint(constraint,target,mesh,fem,filter,mI)
            for i = 1:length(constraint)
                s.type                 = constraint{i};
                s.target               = target;
                s.mesh                 = mesh;
                s.physicalProblem      = fem;
                s.filter               = filter{i};
                s.materialInterpolator = mI;
                sF{i}                  = ShapeFunctional.create(s);
            end
            ss.shapeFunctions = sF;
            sFConstraint      = Constraint(ss);
        end

        function s = createOptimizer(type,primal,monitoring,cost,constraint,x,lam,maxIter,constraintCase,target)
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
            s.primal         = primal;
            switch x.type
                case 'Density'
                    s.ub = 1;
                    s.lb = 0;
            end
            opt = Optimizer.create(s);
            opt.solveProblem();
        end
    end
end