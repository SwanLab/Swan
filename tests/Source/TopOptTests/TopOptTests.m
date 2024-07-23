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
            'test_dirichletProjection','test_gripping','test_micro', 'test_micro2'
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
            bc     = gid.boundaryConditions;
            x      = obj.createDesignVariable(designVariable,m,geomFunSettings,plotting);
            filtersCost = obj.createFilters(filterCostType,m,filterCostSettings);
            filtersConstraint = obj.createFilters(filterConstraintType,m,filterConstraintSettings);
            mI     = obj.createMaterialInterpolator(materialType,method,m,dim);
            mat    = obj.createMaterial(x,mI);
            fem    = obj.createElasticProblem(m,mat,ptype,dim,bc);
            Msmooth = obj.createMassMatrix(m);
            if exist('micro')
                s = micro;
            else
                s = [];
            end
            sFCost = obj.createCost(cost,weights,m,fem,filtersCost,mat,Msmooth,filename,s);
            sFConstraint = obj.createConstraint(constraint,target,m,fem,filtersConstraint,mat,Msmooth);
            l.nConstraints = length(constraint);
            lam    = DualVariable(l);
            primal = optimizerUnconstrained;
            t      = obj.createOptimizer(optimizer,primal,monitoring,sFCost,sFConstraint,x,lam,maxiter,constraint_case,target);
            close all;
        end

        function s = readGidFile(file)
            a.fileName = file;
            s          = FemDataContainer(a);
        end

        function x = createDesignVariable(type,mesh,gSet,plotting)
            g      = GeometricalFunction(gSet);
            lsFun  = g.computeLevelSetFunction(mesh);
            switch type
                case 'Density'
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

        function m = createMaterial(x,mI)
            f = x.obtainDomainFunction();
            f = f.project('P1');
            s.type                 = 'DensityBased';
            s.density              = f;
            s.materialInterpolator = mI;
            s.dim                  = '2D';
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
            s.solverCase         = 'DIRECT';
            s.type               = 'ELASTIC';
            fem                  = PhysicalProblem.create(s);
        end

        function M = createMassMatrix(mesh)
            s.test  = LagrangianFunction.create(mesh,1,'P1');
            s.trial = LagrangianFunction.create(mesh,1,'P1');
            s.mesh  = mesh;
            s.type  = 'MassMatrix';
            LHS = LHSintegrator.create(s);
            M = LHS.compute;
            M = eye(size(M));
        end

        function sFCost = createCost(cost,weights,mesh,fem,filter,mat,Msmooth,filename,s)
            for i = 1:length(cost)
                s.type            = cost{i};
                s.mesh            = mesh;
                s.physicalProblem = fem;
                s.filter          = filter{i};
                s.material        = mat;
                s.filename        = filename;
                sF{i}             = ShapeFunctional.create(s);
            end
            ss.shapeFunctions = sF;
            ss.weights        = weights;
            ss.Msmooth        = Msmooth;
            sFCost            = Cost(ss);
        end

        function sFConstraint = createConstraint(constraint,target,mesh,fem,filter,mat,Msmooth)
            for i = 1:length(constraint)
                s.type            = constraint{i};
                s.target          = target;
                s.mesh            = mesh;
                s.physicalProblem = fem;
                s.filter          = filter{i};
                s.material        = mat;
                s.gradientTest    = LagrangianFunction.create(mesh,1,'P1');
                sF{i}             = ShapeFunctional.create(s);
            end
            ss.shapeFunctions = sF;
            ss.Msmooth        = Msmooth;
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
            s.etaNorm        = 0.05;
            s.gJFlowRatio    = 1; % Only NullSpace
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