classdef VigdergauzParametersFromAxAy < VigdergauzParameters
    
    properties (GetAccess = public, SetAccess = private)
        mx
        my
        FxMax
        FyMax
        R
        rx
        ry
        TxTy
    end
    
    properties (Access = private)
        ax
        ay
        cx
        cy
        Tx
        Ty
    end
       
    methods (Access = public)
        
        function obj = VigdergauzParametersFromAxAy(cParams)
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
            obj.Tx = obj.computeTx();
            obj.Ty = obj.computeTy();
            obj.TxTy = obj.Tx*obj.Ty;
            obj.computeFeasibility()
            obj.rx = obj.computeRx();
            obj.ry = obj.computeRy();            
            obj.R = obj.computeR();           
            obj.FxMax = obj.computeFxMax();
            obj.FyMax = obj.computeFyMax();            
            obj.mx = obj.computeMx();
            obj.my = obj.computeMy();
        end
        
        function Tx = computeTx(obj)
            Tx = obj.computeT(obj.ax,obj.ay,obj.cy);
        end
        
        function Ty = computeTy(obj)
            Ty = obj.computeT(obj.ay,obj.ax,obj.cx);
        end
        
        function computeFeasibility(obj)
            if obj.TxTy > 1 
               error('Tx*Ty > 1'); 
            end
            if obj.TxTy < 0 
               error('Tx*Ty < 0'); 
            end            
        end
        
        function rx = computeRx(obj)
            rx = obj.computeEllipticParameter(obj.Tx);
        end
        
        function ry = computeRy(obj)
            ry = obj.computeEllipticParameter(obj.Ty);            
        end
        
        function x = computeEllipticParameter(obj,T)
            eps = 2.2204*1e-14;
            if T <= 0.15
                x = 1 - eps;
            else
                F = @(x) obj.implicitRequation(x,T);
                x = fzero(F,[0+eps,1-eps]);
            end
        end
              
        function R = computeR(obj)
            f = @(x) (1 - x)/x;
            R = f(obj.rx)*f(obj.ry);
        end    
        
        function Fx = computeFxMax(obj)
            Fx = obj.computeFmax(obj.rx);
        end
        
        function Fy = computeFyMax(obj)
            Fy = obj.computeFmax(obj.ry);
        end
        
        function F = computeFmax(obj,r)
            F = obj.incompleteElliptic(sqrt(1 - obj.R),r);                        
        end
        
        function f = implicitRequation(obj,x,T)
            f = obj.completeElliptic(x)*T - obj.completeElliptic(1-x);
        end
        
        function K = completeElliptic(obj,k)
            K = obj.incompleteElliptic(1,k);
        end
        
        function Mx = computeMx(obj)
            Mx = obj.computeM(obj.ax,obj.rx);            
        end

        function My = computeMy(obj)
            My = obj.computeM(obj.ay,obj.ry);            
        end        
        
        function m = computeM(obj,a,r)
            K = obj.completeElliptic(r);
            Fmax = obj.incompleteElliptic(sqrt(1-obj.R),r);
            m = a/K*Fmax;
        end
        
    end
    
    methods (Access = private, Static)
       
        function F = incompleteElliptic(x,k)
            F = ellipticF(asin(x),k);
        end       
        
        function T = computeT(a1,a2,c)
           T = (c - a2)/a1; 
        end               
        
    end
    
end