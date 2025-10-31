classdef CrossBarsExample < handle

    properties (Access = private)
        mesh
        filter
        designVariable
        materialInterpolator
        physicalProblem
        compliance
        volume
        filterSegment
        perimeter
        cost
        constraint
        primalUpdater
        optimizer
    end

    methods (Access = public)

        function obj = CrossBarsExample()
            obj.init()
            obj.createMesh();
            obj.createDesignVariable();
            obj.createFilter();
            obj.createSegmentFilter();
            obj.createPerimeter();
            obj.perimeter.computeFunctionAndGradient(obj.designVariable)
            obj.designVariable.fun.print('Experiments/LevelSet/CantileverLevelSetGlobalSegmentfValues');
        end

    end

    methods (Access = private)

        function init(obj)
            close all;
        end

        function createMesh(obj)
            obj.mesh = TriangleMesh(1,1,100,50);
        end

        function createDesignVariable(obj)

           
            s.type = 'FourPerpendicularBars';
            s.leftBar_xMax = 0.35;   % right edge of left bar
            s.barWidth = 0.1;

            s.rightBar_xMin = 1 - s.leftBar_xMax;  % left edge of right bar
            s.bottomBar_yMax = s.leftBar_xMax ; % top edge of bottom bar
            s.topBar_yMin = s.rightBar_xMin;    % bottom edge of top bar
            
            
            
            g              = GeometricalFunction(s);
            lsFun          = g.computeLevelSetFunction(obj.mesh);
            s.fun  = lsFun;
            s.mesh = obj.mesh;
            s.type = 'LevelSet';
            s.plotting = true;
            ls     = DesignVariable.create(s);
            obj.designVariable = ls;



        end

        function createFilter(obj)
            s.filterType = 'LUMP';
            s.mesh  = obj.mesh;
            s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            f = Filter.create(s);
            obj.filter = f;
        end

       function uMesh = createBaseDomain(obj)
            levelSet         = -ones(obj.mesh.nnodes,1);
            s.backgroundMesh = obj.mesh;
            s.boundaryMesh   = obj.mesh.createBoundaryMesh();
            uMesh = UnfittedMesh(s);
            uMesh.compute(levelSet);
       end

        function createIsotropicFilter(obj)
            s.mesh  = obj.mesh;
            s.alpha = 4;
            s.beta  = 2*obj.mesh.computeMeanCellSize();
            s.theta = 90;
            s.tol0  = 1e-6;
            obj.filterIso  = NonLinearFilterSegment(s);
        end       

        function createAnisotropicFilter(obj)


        end

        function createSegmentFilter(obj)
            s.mesh  = obj.mesh;
            s.alpha = 4;
            s.beta  = 2*obj.mesh.computeMeanCellSize();
            s.theta = 90;
            s.tol0  = 1e-6;
            obj.filterSegment  = NonLinearFilterSegment(s);
        end

        function createPerimeter(obj)
            h         = obj.mesh.computeMinCellSize();
            s.mesh    = obj.mesh;
            s.uMesh   = obj.createBaseDomain();
            s.filter  = obj.filterSegment;
            s.epsilon = 2*h;
            s.value0  = 1;
            obj.perimeter  = PerimeterFunctional(s);
        end

    end
end