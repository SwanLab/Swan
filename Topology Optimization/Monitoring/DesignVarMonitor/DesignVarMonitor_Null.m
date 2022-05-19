classdef DesignVarMonitor_Null < DesignVarMonitor_Abstract
    
    properties (Access = protected)
        designVarName = []
    end
    
    methods (Access = public)
        
        function obj = DesignVarMonitor_Null(cParams)
            obj@DesignVarMonitor_Abstract(cParams);
            close;
        end
        
        function plot(varargin)
        end
        
        function refresh(obj)
            obj.plot();
            drawnow
        end        
        
    end
    
    methods (Access = protected)
        
        function initPlotting(~)
        end
        
        function createCamera(obj)
            nullAxes = axes;
            obj.cam = Camera_Null(nullAxes);
        end
        
    end
    
    methods (Access = protected, Static)
        
        function color = getColor()
            color = [];
        end
        
    end
    
end