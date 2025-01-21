classdef Camera_Rotatory < Camera_Abstract
    
    properties (Access = protected)
        currentView = [90 30]
    end
    
    properties (Access = private)
        period = 200
        speed
    end
    
    methods (Access = public)
        
        function obj = Camera_Rotatory(axes)
            obj@Camera_Abstract(axes);
        end
        
        function updateView(obj)
            obj.applyAngularMotion();
            obj.checkView();
            obj.refresh();
        end
        
    end
    
    methods (Access = private)
        
        function applyAngularMotion(obj)
            obj.currentView  = obj.currentView + [obj.speed 0];
        end
        
        function checkView(obj)
            az = obj.currentView(1);
            el = obj.currentView(2);
            
            az = obj.keepInRange(az,0,360);
            el = obj.keepInRange(el,-90,90);
            
            obj.currentView(1) = az;
            obj.currentView(2) = el;
            
        end
        
    end
    
    methods (Access = private, Static)
        
        function x = keepInRange(x,xMin,xMax)
            X = abs(xMin-xMax);
            while x > xMax
                x = x - X;
            end
            
            while x < xMin
                x = x + X;
            end
        end
        
    end
    
    methods
        
        function w = get.speed(obj)
            w = 360/obj.period;
        end
        
    end
    
end