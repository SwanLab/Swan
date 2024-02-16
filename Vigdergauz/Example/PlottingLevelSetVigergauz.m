classdef PlottingLevelSetVigergauz < handle
    
    properties (Access = private)
        mesh
        levelSet
        topOpt
    end
    
    methods (Access = public)
        
        function obj = PlottingLevelSetVigergauz()
            obj.init();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.createTopOptProblem();
            obj.createLevelSet();
            obj.createMesh();
            obj.plotMesh();
        end
        
        function createMesh(obj)
            meshBackground = obj.topOpt.designVariable.mesh;
            interpolation = Interpolation.create(meshBackground.type,'LINEAR');
            s.unfittedType = 'INTERIOR';
            s.meshBackground = meshBackground;
            s.interpolationBackground = interpolation;
            s.includeBoxContour = false;
            cParams = SettingsMeshUnfitted(s);
            obj.mesh = UnfittedMesh(cParams);
            obj.mesh.compute(obj.levelSet);
        end
        
        function createLevelSet(obj)
            obj.levelSet = obj.topOpt.designVariable.value;
        end
        
        function createTopOptProblem(obj)
            settings = Settings('VigergauzLevelSetInput');
            translator = SettingsTranslator();
            translator.translate(settings);
            fileName = translator.fileName;
            settingsTopOpt = SettingsTopOptProblem(fileName);
            obj.topOpt = TopOpt_Problem(settingsTopOpt);
        end
               
        function plotMesh(obj)
            obj.mesh.plot();
            view([0 0 1]);
        end
        
    end
    
end