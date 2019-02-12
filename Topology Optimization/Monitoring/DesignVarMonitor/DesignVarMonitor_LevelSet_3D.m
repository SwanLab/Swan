classdef DesignVarMonitor_LevelSet_3D < DesignVarMonitor_LevelSet ...
        %                                       & DesignVarMonitor_3D
    
    properties (Access = protected)
        unfittedType = 'BOUNDARY'
        meshIncludeBoxContour = true
    end
    
    methods (Access = public)
        
        function obj = DesignVarMonitor_LevelSet_3D(mesh)
            obj@DesignVarMonitor_LevelSet(mesh);
        end
        
    end
    
end