classdef testUnfitted < test
    
    properties (Access = protected, Abstract)
        testName
        meshType
        meshIncludeBoxContour
    end
    
    properties (Access = protected)
        topOpt
        levelSet
        unfittedMesh
        oldMeshUnfitted
    end
    
    properties (Access = private)
        settings
    end
    
    methods (Access = protected)
        
        function createTopOpt(obj)
            obj.createSettings();
            obj.topOpt = TopOpt_Problem(obj.settings);
            obj.levelSet = obj.topOpt.designVariable.value;
        end
        
        function createMesh(obj)
            bM  = obj.topOpt.designVariable.mesh;
            s.backgroundMesh = bM.innerMeshOLD;
            s.boundaryMesh   = bM.boxFaceMeshes;
            cParams = SettingsMeshUnfitted(s);            
            obj.unfittedMesh = UnfittedMesh(cParams);
            obj.unfittedMesh.compute(obj.levelSet); 
        end
        
    end
    
    methods (Access = private)
        
        function createSettings(obj)
            obj.createOldSettings();
            obj.translateToNewSettings();
        end
        
        function createOldSettings(obj)
            fileName = obj.testName;
            s = Settings(fileName);
            s.warningHoleBC = false;
            s.printIncrementalIter = false;
            s.printChangingFilter = false;
            s.printing = false;
            obj.settings = s;
        end
        
        function translateToNewSettings(obj)
            translator = SettingsTranslator();
            translator.translate(obj.settings);
            fileName = translator.fileName;
            obj.settings  = SettingsTopOptProblem(fileName);
        end
        
    end
    
end

