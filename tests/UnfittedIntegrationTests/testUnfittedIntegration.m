classdef testUnfittedIntegration < testShowingError
    
    
    properties (Access = protected)
        tol = 6e-2;
        topOpt
        levelSet
    end
    
    properties (Access = protected, Abstract)
        analiticalArea
        meshType
        testName
    end
    
    properties (Access = private)
        areaAdim
    end
    
    methods (Access = protected)
        
        function obj = testUnfittedIntegration()
            obj.createTopOpt()
            obj.integrateSurface()
        end
        
        function createTopOpt(obj)
            file_name_in = strcat('./Input/',obj.testName);
            settings = Settings(file_name_in);
            settings.printChangingFilter = false;
            obj.topOpt = TopOpt_Problem(settings);
            obj.topOpt.preProcess;
        end
        
        function integrateSurface(obj)
            obj.levelSet = obj.topOpt.x;
            obj.createMesh();
            
            area = obj.computeGeometricalVariable();
            obj.areaAdim = area/obj.analiticalArea;
        end
        
        function computeError(obj)
            obj.error = 1 - obj.areaAdim;
        end
    end
    
    methods (Access = protected, Abstract)
        computeGeometricalVariable(obj,levelSet)
        createMesh(obj)
    end
    
end

