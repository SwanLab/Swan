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
            obj.computeLevelSet();
            obj.plotMesh();
        end
        
        function createMesh(obj)
            meshBackground = obj.topOpt.designVariable.mesh;
            interpolation = Interpolation.create(meshBackground,'LINEAR');
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
        
        function computeLevelSet(obj)            
            volume = 0.3;
            theta = 1 - volume;            
            rxMax = 0.99;
            ryMax = 0.99;            
            phimin = atan((theta)/(rxMax^2))*180/pi;
            phimax = atan((rxMax^2)/(theta))*180/pi;
            alpha = 0.2;
            phi = phimin + (phimax - phimin)*alpha;
            phi = phi*pi/180;
            rat = tan(phi);
            x = obj.mesh.meshBackground.coord(:,1);
            y = obj.mesh.meshBackground.coord(:,2);
            obj.levelSet = LevelSetLipung(x,y,theta,phi);
            obj.mesh.compute(obj.levelSet);   
            v = obj.mesh.computeMass();
        end
        
        function plotMesh(obj)
            obj.mesh.plot();
            view([0 0 1]);
        end
        
    end
    
end