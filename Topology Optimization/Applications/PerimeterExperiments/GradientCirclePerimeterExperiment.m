classdef GradientCirclePerimeterExperiment < handle
    
    properties (Access = private)
        radius
        nameCase
        inputFile
        scale
        imesh
        
        mesh
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
                obj.createMesh();
                obj.createLevelSet();
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
        
        function createMesh(obj)
            s.inputFile = obj.inputFile{obj.imesh};
            meshCreator = MeshCreatorFromInputFile(s);
            meshCreator.createMesh();
            obj.mesh         = meshCreator.mesh;
            obj.domainLength = meshCreator.domainLength;
        end
        
        function createLevelSet(obj)
            s.inputFile             = obj.inputFile{obj.imesh};
            s.mesh                  = obj.mesh;
            s.scale                 = obj.scale;
            s.plotting              = true;
            s.printing              = true;
            s.levelSetCreatorParams = obj.createLevelSetCreatorParams();
            lsCreator = LevelSetCreatorForPerimeter(s);
            obj.levelSet = lsCreator.levelSet;
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
            s.mesh             = obj.mesh;
            s.scale            = obj.scale;
            s.designVariable   = obj.levelSet;
            s.outputFigureName = ['SmoothedCircleMesh',num2str(obj.imesh)];
            s.plotting         = false;
            s.printing         = true;
            s.capturingImage   = false;
            s.isRobinTermAdded = false;
            s.perimeterType    = 'perimeterInterior';
            rPerimeter = RegularizedPerimeterComputer(s);
            rPerimeter.compute();
            obj.regularizedPerimeter = rPerimeter;
        end
                
        function computeGradientVariationWithRadius(obj) 
            s.mesh                 = obj.mesh;
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
            s.mesh                 = obj.mesh;
            s.levelSet             = obj.levelSet;
            s.radius               = obj.radius;
            s.regularizedPerimeter = obj.regularizedPerimeter;
            s.domainLength         = obj.domainLength;
            gComputer = GradientVariationWithThetaComputer(s);
            gComputer.compute();
        end        
        
        function computeGradientSurfaces(obj)
            s.outPutFolder = obj.outputFolder;
            s.mesh         = obj.mesh;
            s.rPerimeter   = obj.regularizedPerimeter;
            s.iMesh        = obj.imesh;
            s.domainLength = obj.domainLength;
            gComputer = GradientSurfPerimeterComputer(s);
            gComputer.compute();            
        end
        
    end
    
end