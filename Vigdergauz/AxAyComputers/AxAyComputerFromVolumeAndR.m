classdef AxAyComputerFromVolumeAndR < handle

    properties (Access = private)
        theta
        r
        cx
        cy
    end
    
    methods (Access = public)
        
        function obj = AxAyComputerFromVolumeAndR(cParams)
            obj.init(cParams)
        end
        
        function [ax,ay] = compute(obj)
            h  = obj.computeH(obj.theta,obj.r);
            Tx = obj.computeTx(obj.r,h,obj.cx);
            Ty = obj.computeTy(obj.r,h,obj.cx);
            ax = obj.computeAx(Tx,Ty,obj.cx);
            ay = obj.computeAy(Tx,Ty,obj.cx);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.theta = 1 - cParams.volume;
            obj.r     = cParams.r;
            obj.cx    = cParams.cx;            
            obj.cy    = cParams.cy;
        end
        
    end
    
    methods (Access = private, Static)
        
        function Tx = computeTx(r,h,cx)
            n = r*h*(1+h);
            d = (1+r*h)*cx;
            Tx = n/d;
        end
        
        function Ty = computeTy(r,h,cx)
            n = cx*h*(1+r*h);
            d = (1+h);
            Ty = n/d;
        end
        
        function ax = computeAx(Tx,Ty,cx)
            ax = (cx - Ty)/(1-Tx*Ty);
        end
        
        function ay = computeAy(Tx,Ty,cx)
            ay = (1-Tx*cx)/(1-Tx*Ty);
        end
        
        function h = computeH(theta,r)
            n = -(1+r)*theta + sqrt(theta^2*(r-1)^2 + 4*r);
            d = 2*r*(1+theta);
            h = n/d;
        end
        
    end
    
end
