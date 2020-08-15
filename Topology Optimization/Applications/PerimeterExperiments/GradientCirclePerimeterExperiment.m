classdef GradientCirclePerimeterExperiment < handle
    
    properties (Access = private)
        radius
        nameCase
        inputFile
        scale
        imesh
        
        
        backgroundMesh
        boundaryMesh
        unfittedMesh
        
        domainLength
        levelSet
        
        regularizedPerimeter
        
        xNodesToPlot
        gradientInNodesToPlot
        
        outputFolder

    end
    
    methods (Access = public)
        
        function obj = GradientCirclePerimeterExperiment()
            obj.init();
            for imesh = 1:numel(obj.inputFile)
                obj.imesh = imesh;
                obj.createBackgroundAndBoundaryMesh();
                obj.computeDomainLength();
                obj.createLevelSet();
                obj.createUnfittedMesh();
                obj.computeRegularizedPerimeters();
                obj.computeGradientVariationWithRadius();
                obj.computeGradientVariationWithTheta();
                obj.computeGradientSurfaces();
            end
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.radius = 0.2499999;
            obj.nameCase = 'GradientCirclePerimeterExperiment';
            obj.inputFile = {'SquareMacroTriangle'; ...
                'SquareMacroTriangleFine';...
            'SquareMacroTriangleFineFine'};
            obj.scale     = 'MACRO';
            obj.outputFolder = '/home/alex/Dropbox/Perimeter/';
        end
        
        function createBackgroundAndBoundaryMesh(obj)
            s.inputFile = obj.inputFile{obj.imesh};
            s.isBackgroundMeshRectangularBox = true;
            mCreator = BackgroundAndBoundaryMeshCreatorFromInputFile(s);
            obj.backgroundMesh = mCreator.backgroundMesh;
            obj.boundaryMesh   = mCreator.boundaryMesh;            
        end
        
        
        function computeDomainLength(obj)
            x = obj.backgroundMesh.coord;
            d = max(x(:,1)) - min(x(:,1));          
            obj.domainLength = d;
        end
        
        function createLevelSet(obj)
            s.type = 'circleInclusion';
            halfSide = obj.domainLength/2;            
            s.fracRadius = obj.radius/halfSide;
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
        
        function s = createLevelSetCreatorParams(obj)
            ss.type = 'circleInclusion';
            halfSide = obj.domainLength/2;
            ss.fracRadius = (obj.radius/halfSide);
            s = SettingsLevelSetCreator;
            s = s.create(ss);
        end
        
        function computeRegularizedPerimeters(obj)
            s.inputFile        = obj.inputFile{obj.imesh};
            s.mesh             = obj.backgroundMesh;
            s.scale            = obj.scale;
            s.designVariable   = obj.levelSet;
            s.outputFigureName = ['SmoothedCircleMesh',num2str(obj.imesh)];
            s.plotting         = false;
            s.printing         = false;
            s.capturingImage   = false;
            s.isRobinTermAdded = false;
            s.perimeterType    = 'perimeterInterior';
            rPerimeter = RegularizedPerimeterComputer(s);
            rPerimeter.compute();
            obj.regularizedPerimeter = rPerimeter;
        end
                
        function computeGradientVariationWithRadius(obj) 
            s.mesh                 = obj.backgroundMesh;
            s.regularizedPerimeter = obj.regularizedPerimeter;
            s.inputFile            = obj.inputFile{obj.imesh};
            s.nameCase             = obj.nameCase;
            s.outputFolder         = obj.outputFolder;
            s.domainLength         = obj.domainLength;            
            gComputer = GradientVariationWithRadiusComputer(s);
            gComputer.compute();
        end
        
        function computeGradientVariationWithTheta(obj)
            filePlotName = [obj.outputFolder,obj.inputFile{obj.imesh}];
            s.filePlotName         = filePlotName;
            s.mesh                 = obj.backgroundMesh;
            s.levelSet             = obj.levelSet;
            s.radius               = obj.radius;
            s.regularizedPerimeter = obj.regularizedPerimeter;
            s.domainLength         = obj.domainLength;
            s.circunferenceMesh    = obj.unfittedMesh.boundaryCutMesh;
            gComputer = GradientVariationWithThetaComputer(s);
            gComputer.compute();
        end        
        
        function computeGradientSurfaces(obj)
            s.outPutFolder = obj.outputFolder;
            s.mesh         = obj.backgroundMesh;
            s.rPerimeter   = obj.regularizedPerimeter;
            s.iMesh        = obj.imesh;
            s.domainLength = obj.domainLength;
            gComputer = GradientSurfPerimeterComputer(s);
            gComputer.compute();            
        end
        
    end
    
end