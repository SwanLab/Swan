classdef LevelSetLipung < handle
    
    properties (Access = public)
        value
    end
    
    properties (Access = private)
        x
        y
        volum
        phi
        rx
        ry
        h
        M
        mx
        my
        Tx
        Ty
        r
        ax
        ay
        FxMax
        FyMax
        theta
        cx
        cy
    end
    
    methods (Access = public)
        
        function obj = LevelSetLipung(cParams)
            obj.init(cParams)
            obj.compute();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.x = cParams.x;
            obj.y = cParams.y;
            obj.volum = cParams.volum;
            obj.phi = cParams.phi;
        end
        
        function compute(obj)
            x = obj.x;
            y = obj.y;
            obj.x = (1 - (-1))*(x-0.5);
            obj.y = (1 - (-1))*(y-0.5);
            
            obj.theta = 1-obj.volum;
          
            obj.cx = 1; 
            obj.cy = 1;
            
            obj.r = obj.computeOptimalR();
            
            obj.h = obj.computeH(obj.theta,obj.r);
            
            obj.Tx = obj.computeTx(obj.r,obj.h,obj.cx);
            obj.Ty = obj.computeTy(obj.r,obj.h,obj.cx);
            
            obj.ax = obj.computeAx(obj.Tx,obj.Ty,obj.cx);
            obj.ay = obj.computeAy(obj.Tx,obj.Ty,obj.cx);
            
            obj.mx = obj.computeEllipticParameter(obj.ax,obj.ay,obj.cy);
            obj.my = obj.computeEllipticParameter(obj.ay,obj.ax,obj.cx);
            
            obj.M = obj.computeM(obj.mx,obj.my);
            
            obj.FxMax = obj.incompleteEllipticFunction(sqrt(1 - obj.M),obj.mx);
            obj.FyMax = obj.incompleteEllipticFunction(sqrt(1 - obj.M),obj.my);
            
            obj.rx = obj.computeR(obj.ax,obj.mx,obj.M);
            obj.ry = obj.computeR(obj.ay,obj.my,obj.M);
            
            phi = obj.computeLevelSet(obj.x,obj.y,obj.rx,obj.ry,obj.FxMax,obj.FyMax,obj.mx,obj.my,obj.M);
            obj.value = phi;
        end
        
        function r = computeOptimalR(obj)
            F = @(r) obj.equationForR(r,obj.theta,obj.phi,obj.cx,obj.cy);
            [rub,rlb] = obj.findRbounds2(F);
            r = fzero(F,[rlb,rub]);
        end
        
        function [rub,rlb] = findRbounds2(obj,F)
            r0 = 1;
            F0 = F(1);
            eps = 10^(-6);
            if F0 >= 0
                r1 = 1 + eps;
                F1 = F(r1);
                while F1 >= 0
                    rnew = obj.newPointBySecant(r0,r1,F0,F1);
                    r0 = r1;
                    F0 = F1;
                    r1 = rnew;
                    F1 = F(r1);
                end
                rub = r1;
                rlb = r0;
            else
                r0 = 1;
                r1 = 1 - eps;
                F1 = F(r1);
                while F1 <= 0
                    rnew = obj.newPointBySecant(r0,r1,F0,F1);
                    r0 = r1;
                    F0 = F1;
                    r1 = max(1e-12,rnew);
                    F1 = F(r1);
                end
                rub = r0;
                rlb = r1;
            end
        end
        
        function f = equationForR(obj,r,theta,phi,cx,cy)            
            obj.r = r;
            obj.theta = theta;
            obj.h = obj.computeH(theta,r);
            
            obj.Tx = obj.computeTx(r,obj.h,obj.cx);
            obj.Ty = obj.computeTy(r,obj.h,obj.cx);
            
            obj.ax = obj.computeAx(obj.Tx,obj.Ty,cx);
            obj.ay = obj.computeAy(obj.Tx,obj.Ty,cx);
            
            obj.mx = obj.computeEllipticParameter(obj.ax,obj.ay,obj.cy);
            obj.my = obj.computeEllipticParameter(obj.ay,obj.ax,obj.cx);
            
            obj.M  = obj.computeM(obj.mx,obj.my);
            obj.rx = obj.computeR(obj.ax,obj.mx,obj.M);
            obj.ry = obj.computeR(obj.ay,obj.my,obj.M);
            
            f = tan(phi) - obj.rx/obj.ry;
        end
        
       function r = computeR(obj,a,m,M)
            K = obj.completeEllipticFunction(m);
            Fmax = obj.incompleteEllipticFunction(sqrt(1-M),m);
            r = a/K*Fmax;
        end
        
        function x = computeEllipticParameter(obj,a1,a2,c)
            T = obj.computeT(a1,a2,c);
            eps = 2.2204*1e-14;
            if T <= 0.15
                x = 1 - eps;
            else
                F = @(x) obj.implicitMequation(x,T);
                options = optimset('Display','iter');
                %x = fzero(F,[0+eps,1-eps],options);
                x = fzero(F,[0+eps,1-eps]);
            end
        end
        
        function f = implicitMequation(obj,x,T)
            f = obj.completeEllipticFunction(x)*T - obj.completeEllipticFunction(1-x);
        end
        
        function K = completeEllipticFunction(obj,k)
            K = obj.incompleteEllipticFunction(1,k);
        end
        
    end
    
    methods (Access = private, Static)
        
        function phi = computeLevelSet(x,y,rx,ry,FxMax,FyMax,mx,my,M)
            xp(:,1) = (ellipj(x(:,1)/rx*FxMax,mx));
            yp(:,1) = (ellipj(y(:,1)/ry*FyMax,my));
            levelset(:,1) = (1-xp(:,1).^2).*(1-yp(:,1).^2) - M;
            validx = abs(x) <= rx;
            validy = abs(y) <= ry;
            valid =  validx & validy;
            phi = levelset;
            phi(~(valid)) = -abs(phi(~(valid)));
            
        end
        
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
                
        function x2 = newPointBySecant(x0,x1,f0,f1)
            x2 = x1 - (x1-x0)/(f1 - f0)*f1;
        end
        
        function h = computeH(theta,r)
            n = -(1+r)*theta + sqrt(theta^2*(r-1)^2 + 4*r);
            d = 2*r*(1+theta);
            h = n/d;
        end
        
        function Tx = computeT(a1,a2,c)
            Tx = (c - a2)/a1;
        end
        
        function M = computeM(mx,my)
            f = @(x) (1 - x)/x;
            M = f(mx)*f(my);
        end
        
        function F = incompleteEllipticFunction(x,k)
            F = ellipticF(asin(x),k);
        end        
        
    end
    
end

