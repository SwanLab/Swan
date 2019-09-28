classdef LevelSetLipung < handle
    
    properties (Access = public)
        value
    end
    
    properties (Access = private)
        x
        y
        volum
        phi
        h
        Tx
        Ty
        r
        ax
        ay
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
            obj.x     = cParams.x;
            obj.y     = cParams.y;
            obj.volum = cParams.volum;
            obj.phi   = cParams.phi;
            obj.cx = 1; 
            obj.cy = 1;            
        end
        
        function compute(obj)
            obj.rescaleCoordinates();
            obj.theta = 1-obj.volum;
            obj.r = obj.computeOptimalR();           
            axay = AxAyComputerFromVolumeAndR(obj.volum,obj.r,obj.cx);
            [obj.ax,obj.ay] = axay.compute();
            
            p = computeVigergauzParameters(obj);
            
            obj.value = obj.computeLevelSet(obj.x,obj.y,p);
        end
        
        function rescaleCoordinates(obj)
            obj.x = (1 - (-1))*(obj.x-0.5);
            obj.y = (1 - (-1))*(obj.y-0.5);            
        end
        
        function f = equationForR(obj,r,theta,phi)            
            obj.r = r;
            obj.theta = theta;
            axay = AxAyComputerFromVolumeAndR(obj.volum,r,obj.cx);
            [obj.ax,obj.ay] = axay.compute();
            p = obj.computeVigergauzParameters();
            f = tan(phi) - p.rx/p.ry;
        end
        
        function p = computeVigergauzParameters(obj)
            s.ax = obj.ax;
            s.ay = obj.ay;
            s.cx = obj.cx;
            s.cy = obj.cy;
            p = VigdergauzParametersComputer(s);
        end
        
          function r = computeOptimalR(obj)
            s.x0 = 1;
            s.functionToSolve = @(r) obj.equationForR(r,obj.theta,obj.phi);
            solver = ImplicitEquationSolver(s); 
            r = solver.solve();
        end
        
       function r = computeR(obj,a,m,M)
            K = obj.completeEllipticFunction(m);
            Fmax = obj.incompleteEllipticFunction(sqrt(1-M),m);
            r = a/K*Fmax;
        end

    end
    
    methods (Access = private, Static)
        
        function phi = computeLevelSet(x,y,p)
            xp(:,1) = ellipj(x(:,1)/p.rx*p.FxMax,p.mx);
            yp(:,1) = ellipj(y(:,1)/p.ry*p.FyMax,p.my);
            levelset(:,1) = (1-xp(:,1).^2).*(1-yp(:,1).^2) - p.M;
            validx = abs(x) <= p.rx;
            validy = abs(y) <= p.ry;
            valid =  validx & validy;
            phi = levelset;
            phi(~(valid)) = -abs(phi(~(valid)));           
        end
        

        
    end
    
end

