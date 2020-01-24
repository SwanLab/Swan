classdef ExteriorGradientCirclePerimeterExperiment < handle
    
    properties (Access = private)
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
        
        function obj = ExteriorGradientCirclePerimeterExperiment()
            obj.init();
            for imesh = 1:numel(obj.inputFile)
                obj.imesh = imesh;
                obj.createMesh();
                obj.createLevelSet();
                obj.computeRegularizedPerimeters();
                obj.computeGradientVariationWithTheta();
            end
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.nameCase = 'GradientCirclePerimeterExperiment';
            obj.inputFile = {'CircleMacro'; ...
                'CircleMacroFine';...
                'CircleMacroFineFine'};
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
            s.printing              = false;
            s.levelSetCreatorParams = obj.createLevelSetCreatorParams();
            lsCreator = LevelSetCreatorForPerimeter(s);
            obj.levelSet = lsCreator.levelSet;
        end
        
        function s = createLevelSetCreatorParams(obj)
            ss.type = 'full';
            s = SettingsLevelSetCreator;
            s = s.create(ss);
        end
        
        function computeRegularizedPerimeters(obj)
            s.inputFile        = obj.inputFile{obj.imesh};
            s.mesh             = obj.mesh;
            s.scale            = obj.scale;
            s.designVariable   = obj.levelSet;
            s.outputFigureName = ['SmoothedCircleMesh',num2str(obj.imesh)];
            s.plotting         = true;
            s.printing         = true;
            s.capturingImage   = false;
            s.isRobinTermAdded = true;
            s.perimeterType    = 'perimeter';
            rPerimeter = RegularizedPerimeterComputer(s);
            rPerimeter.compute();
            obj.regularizedPerimeter = rPerimeter;
        end
                
        function computeGradientVariationWithTheta(obj)
            filePlotName = [obj.outputFolder,obj.inputFile{obj.imesh}];
            s.filePlotName         = filePlotName;
            s.mesh                 = obj.mesh;
            s.levelSet             = obj.levelSet;
            s.radius               = 1;
            s.regularizedPerimeter = obj.regularizedPerimeter;
            s.domainLength         = obj.domainLength;
            gComputer = ExteriorGradientVariationWithThetaComputer(s);
            gComputer.compute();
        end        
        
    end
    
end