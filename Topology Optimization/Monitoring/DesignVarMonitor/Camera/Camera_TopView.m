classdef Camera_TopView < Camera_Abstract
    
    properties (Access = protected)
        currentView = [0 90]
    end
    
    methods (Access = public)
        
        function obj = Camera_TopView(axes)
           obj@Camera_Abstract(axes);
        end
        
        function updateView(obj)
            obj.refresh();
        end
        
    end
    
end