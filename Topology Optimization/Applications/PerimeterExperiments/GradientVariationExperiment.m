classdef GradientVariationExperiment < handle
    
    properties (Access = public)
        backgroundMesh
        boundaryMesh
        levelSet
        regularizedPerimeter
        domainLength
    end
    
    properties (Access = private)
        iMesh        
        inputFile
        levelSetParams
    end    
  
    methods (Access = public)
        
        function obj = GradientVariationExperiment(cParams)
            obj.init(cParams)
            obj.createBackgroundAndBoundaryMesh();
            obj.computeDomainLength();
            obj.createLevelSet();
            obj.createRegularizedPerimeters();
        end
        
    end
    
    methods (Access = protected)
        
        function init(obj,cParams)
            obj.inputFile      = cParams.inputFile;
            obj.iMesh          = cParams.iMesh;
            obj.levelSetParams = cParams.levelSetParams;
        end
        
        function createBackgroundAndBoundaryMesh(obj)
            s.inputFile = obj.inputFile;
            s.isBackgroundMeshRectangularBox = true;
            mCreator = BackgroundAndBoundaryMeshCreatorFromInputFile(s);
            obj.backgroundMesh = mCreator.backgroundMesh;
            obj.boundaryMesh = mCreator.boundaryMesh;
        end
        
        function createLevelSet(obj)
            s = obj.levelSetParams;
            s.coord      = obj.backgroundMesh.coord;
            s.ndim       = obj.backgroundMesh.ndim;
            lsCreator = LevelSetCreator.create(s);
            obj.levelSet = lsCreator.getValue();
        end
        
        function computeDomainLength(obj)
            x = obj.backgroundMesh.coord;
            d = max(x(:,1)) - min(x(:,1));
            obj.domainLength = d;
        end
        
        function createRegularizedPerimeters(obj)
            s.inputFile        = obj.inputFile;
            s.backgroundMesh   = obj.backgroundMesh;
            s.scale            = 'MACRO';
            s.designVariable   = obj.levelSet;
            s.outputFigureName = ['SmoothedCircleMesh',num2str(obj.iMesh)];
            s.plotting         = false;
            s.printing         = false;
            s.capturingImage   = false;
            s.isRobinTermAdded = true;
            s.perimeterType    = 'perimeterInterior';
            rPerimeter = RegularizedPerimeterComputer(s);
            rPerimeter.compute();
            obj.regularizedPerimeter = rPerimeter;
        end
        
    end
    
end