classdef SolutionSaverInMatFormatFromResFile < handle
    
    properties (Access = private)
        path
        testName
        finalIter
    end
    
    methods (Access = public)
        
        function obj = SolutionSaverInMatFormatFromResFile()
            obj.init()
            obj.save();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.path = '/home/alex/git-repos/Swan/Output/';
            obj.testName = 'CantileverFlat';
        end
        
        function save(obj)
            obj.saveDesignVariable()
            obj.saveUnfittedMesh();
        end
        
        function saveDesignVariable(obj)
            s.testPath  = fullfile(obj.path,obj.testName);
            s.testName  = obj.testName;
            d = DesignVariableAndMeshSaverFromRes(s);
            d.save();
        end
        
        function saveUnfittedMesh(obj)
            s.testPath  = fullfile(obj.path,obj.testName);
            d = UnfittedMeshSaverFromSavedMeshAndLevelSet(s);
            d.save();
        end
        
    end
    
end