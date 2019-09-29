classdef VigdergauzParametersComputerFromAxAy < handle
    
    properties (GetAccess = public, SetAccess = private)
        mx
        my
        FxMax
        FyMax
        M
        rx
        ry
    end
    
    properties (Access = private)
        ax
        ay
        cx
        cy
    end
       
    methods (Access = public)
        
        function obj = VigdergauzParametersComputerFromAxAy(cParams)
            obj.init(cParams)
            obj.computeParameters();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.ax = cParams.ax;
            obj.ay = cParams.ay;
            obj.cx = cParams.cx;
            obj.cy = cParams.cy;
        end
        
        function computeParameters(obj)
            obj.mx = obj.computeMx();
            obj.my = obj.computeMy();            
            obj.M = obj.computeM();           
            obj.FxMax = obj.computeFxMax();
            obj.FyMax = obj.computeFyMax();            
            obj.rx = obj.computeRx();
            obj.ry = obj.computeRy();
        end
        
        function mx = computeMx(obj)
            mx = obj.computeEllipticParameter(obj.ax,obj.ay,obj.cy);
        end
        
        function my = computeMy(obj)
            my = obj.computeEllipticParameter(obj.ay,obj.ax,obj.cx);            
        end
        
        function x = computeEllipticParameter(obj,a1,a2,c)
            T = (c - a2)/a1;
            eps = 2.2204*1e-14;
            if T <= 0.15
                x = 1 - eps;
            else
                F = @(x) obj.implicitMequation(x,T);
                x = fzero(F,[0+eps,1-eps]);
            end
        end
              
        function M = computeM(obj)
            f = @(x) (1 - x)/x;
            M = f(obj.mx)*f(obj.my);
        end    
        
        function Fx = computeFxMax(obj)
            Fx = obj.computeFmax(obj.mx);
        end
        
        function Fy = computeFyMax(obj)
            Fy = obj.computeFmax(obj.my);
        end
        
        function F = computeFmax(obj,m)
            F = obj.incompleteElliptic(sqrt(1 - obj.M),m);                        
        end
        
        function f = implicitMequation(obj,x,T)
            f = obj.completeElliptic(x)*T - obj.completeElliptic(1-x);
        end
        
        function K = completeElliptic(obj,k)
            K = obj.incompleteElliptic(1,k);
        end
        
        function Rx = computeRx(obj)
            Rx = obj.computeR(obj.ax,obj.mx);            
        end

        function Ry = computeRy(obj)
            Ry = obj.computeR(obj.ay,obj.my);            
        end        
        
        function r = computeR(obj,a,m)
            K = obj.completeElliptic(m);
            Fmax = obj.incompleteElliptic(sqrt(1-obj.M),m);
            r = a/K*Fmax;
        end
        
    end
    
    methods (Access = private, Static)
       
        function F = incompleteElliptic(x,k)
            F = ellipticF(asin(x),k);
        end                 
        
    end
    
end