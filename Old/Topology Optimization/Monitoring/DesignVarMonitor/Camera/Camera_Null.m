classdef Camera_Null < Camera_Abstract
    
    properties (Access = protected)
        currentView = []
    end
    
    methods (Access = public)
        
        function obj = Camera_Null(axes)
           obj@Camera_Abstract(axes);
        end
        
        function updateView(obj)
        end
        
    end
    
end