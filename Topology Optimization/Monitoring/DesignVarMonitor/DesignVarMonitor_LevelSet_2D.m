classdef DesignVarMonitor_LevelSet_2D < DesignVarMonitor_LevelSet ...
        %                                       & DesignVarMonitor_2D
    
    properties (Access = protected)
        unfittedType = 'INTERIOR'
        meshIncludeBoxContour = false
    end
    
    methods (Access = public)
        
        function obj = DesignVarMonitor_LevelSet_2D(mesh)
            obj@DesignVarMonitor_LevelSet(mesh);
        end
        
    end
    
end