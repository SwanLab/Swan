classdef LocalPerimeterValidation < handle

    properties (Access = private)
        mesh
        designVariable
        perimeter
        volume
        cost
        constraint
        primalUpdater
        optimizer
    end

    methods (Access = public)

        function obj = LocalPerimeterValidation()
            obj.init()
            obj.createMesh();
            obj.createDesignVariable();
            obj.createPerimeter();
            obj.createVolumeConstraint();
            obj.createCost();
            obj.createConstraint();
            obj.createPrimalUpdater();
            obj.createOptimizer();

            saveas(gcf,'Paper/Validation/MonitoringLocalPerimeterValidation.fig');
            obj.designVariable.fun.print('Paper/Validation/LocalPerimeterValidationfValues');
        end

    end

    methods (Access = private)

        function init(obj)
            close all;
        end

        function createMesh(obj)
            obj.mesh = TriangleMesh(3,3,50,50);
        end

        function createDesignVariable(obj)
            s.type  = 'Square';
            s.length = 0.4;
            xC = [0.5,1.5,2.5,0.5,1.5,2.5,0.5,1.5,2.5];
            yC = [2.5,2.5,2.5,1.5,1.5,1.5,0.5,0.5,0.5];

            lsMat = zeros(obj.mesh.nnodes,length(xC));
            for i = 1:length(xC)
                s.xCoorCenter = xC(i);
                s.yCoorCenter = yC(i);
                g             = GeometricalFunction(s);
                ls            = g.computeLevelSetFunction(obj.mesh);
                lsMat(:,i)    = ls.fValues;
            end

            lsVal = min(lsMat,[],2);
            lsFun = LagrangianFunction.create(obj.mesh,1,'P1');
            lsFun.setFValues(-lsVal);

            s.fun              = lsFun;
            s.mesh             = obj.mesh;
            s.type             = 'LevelSet';
            s.plotting         = true;
            ls                 = DesignVariable.create(s);
            obj.designVariable = ls;
        end

        function uMesh = createBaseDomainVolumeCost(obj)
            levelSet         = -ones(obj.mesh.nnodes,1);
            s.backgroundMesh = obj.mesh;
            s.boundaryMesh   = obj.mesh.createBoundaryMesh();
            uMesh = UnfittedMesh(s);
            uMesh.compute(levelSet);
        end

        function uMesh = createBaseDomainPerimeterConstraint(obj,x0,y0)
            s.type             = 'Square';
            s.xCoorCenter      = x0;
            s.yCoorCenter      = y0;
            s.length           = 1;
            g                  = GeometricalFunction(s);
            lsFun              = g.computeLevelSetFunction(obj.mesh);
            sUm.backgroundMesh = obj.mesh;
            sUm.boundaryMesh   = obj.mesh.createBoundaryMesh();
            uMesh              = UnfittedMesh(sUm);
            uMesh.compute(lsFun.fValues);
        end

        function createPerimeter(obj)
            sF.mesh       = obj.mesh;
            sF.filterType = 'PDE';
            sF.trial      = LagrangianFunction.create(obj.mesh,1,'P1');
            fPDE          = Filter.create(sF);
            f{1} = fPDE;
            f{2} = fPDE;
            f{3} = fPDE;

            CAnisotropic = [tand(80), 0; 0, 1/tand(80)];
            aniAlphaDeg = 90;
            R = [cosd(aniAlphaDeg),-sind(aniAlphaDeg)
                sind(aniAlphaDeg), cosd(aniAlphaDeg)];
            CGlobal = R*CAnisotropic*R';
            sF.boundaryType = 'Neumann';
            sF.metric       = 'Anisotropy';
            sF.trial        = LagrangianFunction.create(obj.mesh,1,'P1');
            sF.A            = ConstantFunction.create(CGlobal,obj.mesh);
            fAni1            = Filter.create(sF);
            f{4} = fAni1;

            sF.alpha = 4;
            sF.beta  = 0;
            sF.theta = 90;
            sF.tol0  = 1e-6;
            fSeg    = NonLinearFilterSegment(sF);
            f{5} = fSeg;

            aniAlphaDeg = 90;
            R = [cosd(aniAlphaDeg),-sind(aniAlphaDeg)
                sind(aniAlphaDeg), cosd(aniAlphaDeg)];
            CGlobal = R*CAnisotropic*R';
            sF.trial        = LagrangianFunction.create(obj.mesh,1,'P1');
            sF.A            = ConstantFunction.create(CGlobal,obj.mesh);
            fAni2            = Filter.create(sF);
            f{6} = fAni2;

            aniAlphaDeg = 90+45;
            R = [cosd(aniAlphaDeg),-sind(aniAlphaDeg)
                sind(aniAlphaDeg), cosd(aniAlphaDeg)];
            CGlobal = R*CAnisotropic*R';
            sF.trial        = LagrangianFunction.create(obj.mesh,1,'P1');
            sF.A            = ConstantFunction.create(CGlobal,obj.mesh);
            fAni3            = Filter.create(sF);
            f{7} = fAni3;

            aniAlphaDeg = 0;
            R = [cosd(aniAlphaDeg),-sind(aniAlphaDeg)
                sind(aniAlphaDeg), cosd(aniAlphaDeg)];
            CGlobal = R*CAnisotropic*R';
            sF.trial        = LagrangianFunction.create(obj.mesh,1,'P1');
            sF.A            = ConstantFunction.create(CGlobal,obj.mesh);
            fAni4            = Filter.create(sF);
            f{8} = fAni4;

            aniAlphaDeg = 45;
            R = [cosd(aniAlphaDeg),-sind(aniAlphaDeg)
                sind(aniAlphaDeg), cosd(aniAlphaDeg)];
            CGlobal = R*CAnisotropic*R';
            sF.trial        = LagrangianFunction.create(obj.mesh,1,'P1');
            sF.A            = ConstantFunction.create(CGlobal,obj.mesh);
            fAni5            = Filter.create(sF);
            f{9} = fAni5;

            h         = obj.mesh.computeMeanCellSize();
            s.mesh    = obj.mesh;
            s.epsilon = 2*h;
            s.minEpsilon = 2*h;
            s.value0 = 1;
            s.tarVolume = 0.4;

            tarRef = [0.65,1,1.35, 0.65,1,1.15, 1,1,1];
            x0 = repmat([0.5,1.5,2.5],[1,3]);
            y0 = [repmat(2.5,[1,3]),repmat(1.5,[1,3]),repmat(0.5,[1,3])];
            for i = 1:length(x0)
                s.uMesh          = obj.createBaseDomainPerimeterConstraint(x0(i),y0(i));
                s.target         = tarRef(i);
                s.target0        = s.target;
                s.filter         = f{i};
                obj.perimeter{i} = PerimeterConstraint(s);
            end
        end

        function createVolumeConstraint(obj)
            s.mesh   = obj.mesh;
            s.test = LagrangianFunction.create(obj.mesh,1,'P1');
            s.uMesh = obj.createBaseDomainVolumeCost();
            v = VolumeFunctional(s);
            obj.volume = v;
        end

        function createCost(obj)
            s.shapeFunctions{1} = obj.volume;
            s.weights           = 1;
            s.Msmooth           = obj.createMassMatrix();
            obj.cost            = Cost(s);
        end

        function M = createMassMatrix(obj)
            n = obj.mesh.nnodes;
            h = obj.mesh.computeMinCellSize();
            M = h^2*sparse(1:n,1:n,ones(1,n),n,n);
        end

        function createConstraint(obj)
            for i = 1:length(obj.perimeter)
                s.shapeFunctions{i} = obj.perimeter{i};
            end
            s.Msmooth           = obj.createMassMatrix();
            obj.constraint      = Constraint(s);
        end

        function createPrimalUpdater(obj)
            s.mesh = obj.mesh;
            obj.primalUpdater = SLERP(s);
        end

        function createOptimizer(obj)
            s.monitoring     = true;
            s.cost           = obj.cost;
            s.constraint     = obj.constraint;
            s.designVariable = obj.designVariable;
            s.maxIter        = 100;
            s.tolerance      = 1e-8;
            s.constraintCase = repmat({'EQUALITY'},[1,9]);
            s.etaNorm        = 0.01;
            s.etaNormMin     = 0.01;
            s.gJFlowRatio    = 8.0;
            s.etaMax         = 100;
            s.etaMaxMin      = 0.1;
            s.primalUpdater  = obj.primalUpdater;
            s.gif            = false;
            s.gifName        = [];
            s.printing       = false;
            s.printName      = [];
            opt = OptimizerNullSpace(s);
            opt.solveProblem();
            obj.optimizer = opt;
        end
    end
end