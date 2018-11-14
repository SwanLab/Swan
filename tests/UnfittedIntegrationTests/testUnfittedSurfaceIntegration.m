classdef testUnfittedSurfaceIntegration < testShowingError ...
            & testLoadStoredVariable ...
            & testStoredComputedChecker
   
        
   properties (Access = protected)
      variablesToStore = {'A_star_ref'};
      tol = 1e-9;      
   end
        
   properties (Access = private)
       topOpt
       areaAdim
   end
    
   methods (Access = protected)
       
       function obj = testUnfittedSurfaceIntegration()
           obj.createTopOpt()
           obj.integrateSurface()
           obj.computeAnalyticalArea()
       end
       
        function createTopOpt(obj)
            file_name_in = strcat('./Input/',obj.testName);
            settings = Settings(file_name_in);
            settings.printChangingFilter = false;
            obj.topOpt = TopOpt_Problem(settings);
            obj.topOpt.preProcess;
        end
        
        function integrateSurface(obj)
            settings = obj.topOpt.settings;
            filter =  Filter.create(settings);
            filter.preProcess;
            levelSet = obj.topOpt.x;
            area = filter.computeFacetSurface(levelSet);
            analiticalArea = 4*pi;
            obj.areaAdim = area/analiticalArea;            
        end
        
        function computeAnalyticalArea(obj)
            obj.computedVar{1} = obj.areaAdim;            
        end
        
    end
    
end

