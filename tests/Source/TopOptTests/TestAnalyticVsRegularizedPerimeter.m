classdef TestAnalyticVsRegularizedPerimeter < handle
    
    properties (Access = private)
        designVariable
        backgroundMesh
        boundaryMesh
     %   unfittedMesh
        levelSet
        regularizedPerimeter
        analyticPerimeter
        radius
        inputFile
    end
    
    properties (Access = protected)
        testName = 'testAnalyticVsRegularizedPerimeter';
    end
    
    methods (Access = public)
        
        function obj = TestAnalyticVsRegularizedPerimeter()
            obj.init();
            obj.createBackgroundAndBoundaryMesh();
            obj.createLevelSet();
            obj.computeRegularizedPerimeter();
        end

        function error = computeError(obj)
            aP = obj.analyticPerimeter;
            nP = obj.regularizedPerimeter;
            error = abs(aP - nP)/aP;
        end

    end

    methods (Access = private)
        
        function init(obj)
            obj.radius = 0.31;
            obj.analyticPerimeter = 2*pi*obj.radius;
            obj.inputFile = 'SquareMacroTriangle';
        end
        
        function createBackgroundAndBoundaryMesh(obj)
            s.inputFile = obj.inputFile;
            s.isBackgroundMeshRectangularBox = true;
            mCreator = BackgroundAndBoundaryMeshCreatorFromInputFile(s);
            obj.backgroundMesh = mCreator.backgroundMesh;
            obj.boundaryMesh   = mCreator.boundaryMesh;
        end
        
        function createLevelSet(obj)
            domainLength = max(obj.backgroundMesh.coord(:,1)) - min(obj.backgroundMesh.coord(:,1));
            halfSide = domainLength/2;
            sa.fracRadius = (obj.radius/halfSide);
            sa.coord      = obj.backgroundMesh.coord;
            sa.ndim       = obj.backgroundMesh.ndim;
            s.creatorSettings = sa;
            s.initialCase = 'circleInclusion';
            s.type = 'LevelSet';
            s.mesh = Mesh_Total(obj.backgroundMesh);
            s.scalarProductSettings.epsilon = 0.01;
            s.scalarProductSettings.mesh = Mesh_Total(obj.backgroundMesh);
            s.isFixed = [];
            s.value = [];
            obj.levelSet = DesignVariable.create(s);
        end

        function computeRegularizedPerimeter(obj)
            s = obj.createPerimeterParams();
            perimeterFunc = ShFunc_Perimeter(s);
            perimeterFunc.computeFunctionAndGradient();
            per = perimeterFunc.value*perimeterFunc.value0;
            obj.regularizedPerimeter = per;
        end
        
        function s = createPerimeterParams(obj)
           sC.inputFile      = obj.inputFile;
           sC.mesh           = obj.backgroundMesh;
           sC.designVariable = obj.levelSet;
           sC.epsilon        = obj.backgroundMesh.computeMeanCellSize();
           sC.scale          = 'MACRO';
           sC.type           = 'perimeterInterior';
           sC.isRobinTermAdded = false;
           fCreator = PerimeterParamsCreator(sC);
           s = fCreator.perimeterParams;
        end

    end

end