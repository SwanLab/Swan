classdef DesignVarMonitor_LevelSet_2D < DesignVarMonitor_LevelSet
    
    properties (Access = protected)
        unfittedType = 'INTERIOR'
        meshIncludeBoxContour = false
    end
    
    methods (Access = public)
        
        function obj = DesignVarMonitor_LevelSet_2D(deisgnVar,showBC)
            obj@DesignVarMonitor_LevelSet(deisgnVar,showBC);
        end
        
    end
    
end