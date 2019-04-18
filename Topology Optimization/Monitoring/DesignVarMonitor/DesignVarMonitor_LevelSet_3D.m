classdef DesignVarMonitor_LevelSet_3D < DesignVarMonitor_LevelSet
    
    properties (Access = protected)
        unfittedType = 'BOUNDARY'
        meshIncludeBoxContour = true
    end
    
    methods (Access = public)
        
        function obj = DesignVarMonitor_LevelSet_3D(designVar,showBC)
            obj@DesignVarMonitor_LevelSet(designVar,showBC);
        end
        
    end
    
end