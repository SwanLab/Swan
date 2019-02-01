classdef DesignVarMonitor_Density_3D < DesignVarMonitor_Density %...
                                       % & DesignVarMonitor_3D
    
    methods (Access = public)
        
        function obj = DesignVarMonitor_Density_3D(mesh)
            obj@DesignVarMonitor_Density(mesh);
        end
        
    end
    
    methods (Access = protected)
        
        function init(obj)
            
        end
        
    end
    
end