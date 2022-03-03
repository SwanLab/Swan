classdef UnitBallProjector < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = UnitBallProjector(cParams)
            obj.init(cParams)            
        end
        
        function vN = project(obj,v)
            vx = v(:,1);
            vy = v(:,2);
            norm = sqrt(vx.^2 + vy.^2);
            vN(:,1) = vx./norm;
            vN(:,2) = vy./norm;            
        end
        
        function optD = computeDualOptimality(obj,v)
            vx = v(:,1);
            vy = v(:,2);
            normV = sqrt(vx.^2 + vy.^2);
            optD = norm(normV - 1);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            
        end
        
    end
    
end