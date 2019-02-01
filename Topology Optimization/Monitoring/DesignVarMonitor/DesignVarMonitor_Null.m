classdef DesignVarMonitor_Null < DesignVarMonitor_Abstract
    
    properties (Access = protected)
        designVarible = []
    end
    
    methods (Access = public)
        
        function obj = DesignVarMonitor_Null(mesh)
            obj@DesignVarMonitor_Abstract(mesh);
            close;
        end
        
        function plot(~)
        end
        
    end
    
    methods (Access = protected)
        
        function init(~)
        end
        
    end
    
    methods (Access = protected, Static)
        
        function color = getColor()
            color = [];
        end
        
    end
    
end