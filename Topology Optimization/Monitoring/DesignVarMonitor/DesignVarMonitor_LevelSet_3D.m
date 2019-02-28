classdef DesignVarMonitor_LevelSet_3D < DesignVarMonitor_LevelSet
    
    properties (Access = protected)
        unfittedType = 'BOUNDARY'
        meshIncludeBoxContour = true
    end
    
    methods (Access = public)
        
        function obj = DesignVarMonitor_LevelSet_3D(mesh,showBC)
            obj@DesignVarMonitor_LevelSet(mesh,showBC);
        end
        
    end
    
end