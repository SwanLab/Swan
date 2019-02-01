classdef DesignVarMonitor_Density < DesignVarMonitor_Abstract
    
    properties (Access = protected)
        designVarName = 'Density - \rho';
    end
    
    methods (Access = public)
        
        function obj = DesignVarMonitor_Density(mesh)
            obj@DesignVarMonitor_Abstract(mesh);
        end
        
        function plot(obj,rho)
            set(obj.patchHandle,'FaceVertexAlphaData',double(rho));
            view([1 1 1])
        end
        
    end
    
    methods (Access = protected,Static)
        
        function color = getColor()
            color = [0 0 0];
        end
        
    end
    
end