classdef EdcByLevelSet < ElementalDensityCreator
    
    properties (Access = private)
       levelSet
    end
    
    methods (Access = public)
       
        function ls = getLevelSet(obj)
            ls = obj.levelSet;
        end
    end
    
    methods (Access = protected)
        
        function createDensity(obj,d)
            obj.createLevelSet(d);
            obj.computeDensity(d);
        end
        
    end
    
    methods (Access = private)
       
        function createLevelSet(obj,d)
           obj.levelSet = d.levelSet;
        end
        
        function computeDensity(obj,d)
            lS = obj.levelSet;
            d = d.filterDataBase;
            filter = FilterP0(lS,d);
            obj.density = filter.getDensity();
        end
        
    end
    
end