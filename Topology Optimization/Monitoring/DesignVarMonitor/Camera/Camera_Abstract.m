classdef Camera_Abstract < handle
    
    properties(Access = protected, Abstract)
        currentView
    end
    
    properties (Access = protected)
        axes
    end
    
    methods (Access = public, Abstract)
        
        updateView(obj)
        
    end
    
    methods (Access = public)
        
        function obj = Camera_Abstract(axes)
            obj.axes = axes;
        end
        
    end
    
    methods (Access = protected)
        
        function obj = refresh(obj)
            obj.axes.View = obj.currentView;
        end
        
    end
    
end