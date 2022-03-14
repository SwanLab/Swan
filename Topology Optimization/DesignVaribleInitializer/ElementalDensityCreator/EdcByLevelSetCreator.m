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
        
        function createDensity(obj,cParams)
            obj.createLevelSet(cParams);
            obj.computeDensity(cParams);
        end
        
    end
    
    methods (Access = private)
       
        function createLevelSet(obj,cParams)
           d = cParams.levelSetCreatorDataBase;
           lsC = LevelSetCreator.create(d);
           obj.levelSet = lsC.getValue();
        end
        
        function computeDensity(obj,cParams)
            lS = obj.levelSet;
            d = cParams.filterDataBase;
            filter = FilterP0(lS,d);
            obj.density = filter.getDensity();
        end
        
    end
    
    
end