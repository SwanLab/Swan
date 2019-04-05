classdef testExplorSettLevSetCustom < testShowingError & ...
        testLoadStoredVariable & ...
        testStoredComputedChecker
    
    properties (Access = protected)
        tol = 1e-13;
        testName = 'testExplorSettLevSetCustom';
        variablesToStore = {'levelSet'};
    end
    
    properties (Access = private)
        levelSet
        levelSetParams
    end
    
    methods (Access = public)
        
        function obj = testExplorSettLevSetCustom()
            obj.init();
            obj.computeLevelSet();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.levelSetParams = SettingsLevelSetCreator('paramsLevelSetCreator_circle');
            coord = Mesh_GiD('test2d_quad.m').coord;
            obj.levelSetParams.coord = coord;
        end
        
        function computeLevelSet(obj)
            d = obj.levelSetParams;
            lsC = LevelSetCreator.create(d);
            obj.levelSet = lsC.getValue();
            obj.computedVar{1} = obj.levelSet;
        end
        
    end
    
end