classdef testAnalyticVsRegularizedPerimeter < testShowingError
    
    properties (Access = private)
        designVariable
        backgroundMesh
        boundaryMesh
        unfittedMesh
        levelSet
        regularizedPerimeter
        analyticPerimeter
        radius
        inputFile
    end
    
    properties (Access = protected)
        testName = 'testAnalyticVsRegularizedPerimeter';
        tol
    end
    
    methods (Access = public)
        
        function obj = testAnalyticVsRegularizedPerimeter()
            obj.init();
            obj.createBackgroundAndBoundaryMesh();            
            obj.createLevelSet();
            obj.createUnfittedMesh();
            obj.computeRegularizedPerimeter();
        end
        
    end
    
    methods (Access = protected)
        
        function computeError(obj)
            aP = obj.analyticPerimeter;
            nP = obj.regularizedPerimeter;
            obj.error = abs(aP - nP)/aP;
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.tol = 5e-2;
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
            s.type = 'circleInclusion';
            
            domainLength = max(obj.backgroundMesh.coord(:,1)) - min(obj.backgroundMesh.coord(:,1));
            halfSide = domainLength/2;
            
            s.fracRadius = (obj.radius/halfSide);
            s.coord      = obj.backgroundMesh.coord;
            s.ndim       = obj.backgroundMesh.ndim;
            lsCreator = LevelSetCreator.create(s);
            obj.levelSet = lsCreator.getValue();
        end        
        
        function createUnfittedMesh(obj)
            s.backgroundMesh = obj.backgroundMesh;
            s.boundaryMesh   = obj.boundaryMesh;
            uMesh = UnfittedMesh(s);
            uMesh.compute(obj.levelSet);  
            obj.unfittedMesh = uMesh;
        end        
          
         function computeRegularizedPerimeter(obj)
            s = obj.createPerimeterParams();
            perimeterFunc = ShFunc_Perimeter(s);
            perimeterFunc.computeCostAndGradient();
            obj.regularizedPerimeter = perimeterFunc.value;             
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