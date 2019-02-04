classdef EdcByLevelSetCreator < ElementalDensityCreator
    
    properties (Access = private)
       levelSet
    end
    
    methods (Access = public)
       
        function ls = getLevelSet(obj)
            ls = obj.levelSet;
        end
        
        function f = getFieldsToPrint(obj)
            f{1} = obj.density();
            f{2} = obj.levelSet();
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
           d = d.levelSetCreatorDataBase;             
           lsC = LevelSetCreator.create(d);
           obj.levelSet = lsC.getValue();             
        end
        
        function computeDensity(obj,d)
            lS = obj.levelSet;
            d = d.filterDataBase; 
            filter = FilterP0(lS,d);
            obj.density = filter.getDensity();            
        end
        
    end
    
    
end