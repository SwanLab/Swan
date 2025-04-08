classdef Kernel
    % Kernel
    % Represents the kernel G(x)
    properties (SetAccess = protected, GetAccess = public)
        func, % anonymous function representing G(x)
        der, % anonymous function representing G'(x)
        scalFunc, % anonymous function, 
        % scalFunc(a,b,rho) =2 \pi C(rho) \int_{a}^b rG'(r) J_0(\rho r)
        % where C(rho) is the normalization constant, si Cp.m in folder
        % radialQuad.
        normFunc,  % anonymous function, 
        % normFunc(a,b) = \sqrt{2\pi \int_{a}^b r G'(r)^2 dr}
        startFreq = 0, % real number 
        % frequency where the greatest Bessel coefficients are located. 
        gamma_est, % Estimator used to predict the number of coefficients in 
        % the Bessel series
        customRadialQuad = [],
        % If non-empty, contains the method to compute the radial
        % quadrature optimized for this kernel. 
        singular = true;
        % Set to false if x G'(x)^2 is integrable near 0. 
        lim0;
    end
    
    methods
        % Class constructor
        function[kernel] = Kernel(func,der)
            if nargin == 0
                func = @(x)(0*x);
                der = @(x)(0*x);
            end                
            kernel.func = func;
            kernel.der = der;
            fun1 = @(rho)(@(x)((x.*der(x)).*...
                -rho(:)'.*Cp(rho(:)').*besselj(1,rho(:)'*x)));
            % H10 scalar product
            fun2 = @(x)(x.*der(x).^2);
            % H10 norm
            kernel.scalFunc = @(a,b,rho)(2*pi*integral(fun1(rho),...
                a,b,'ArrayValued',true)'); % -2\pi \int_{a}^b rf'(r)ep'(r)dr
            kernel.normFunc = @(a,b)(sqrt(2*pi*integral(fun2,a,b)));
            kernel.startFreq = 0; % We don't know a priori where the energy
            % is located in the spectrum. 
            kernel.gamma_est = @(tol)deal(0,7); % No fine tuning of gamma
            kernel.customRadialQuad = [];
            if ~ or(isinf(func(0)),isnan(func(0)))
                kernel.singular = false;
                kernel.gamma_est = @(tol)deal(0,3); % No fine tuning of gamma
                kernel.lim0 = func(0);
            else
                kernel.lim0 = 0;
            end
        end
    end
    
    methods (Access = public)
        
        % Efficient Bessel Decomposition for the kernel.        
        function[onlineEBD,rq,loc] = offlineEBD(this,X,Y,a,tol)
            
            % Rescaling
            rMax = rMaxCalc(X,Y); % diameter
            x = X/rMax; % in B(0,1)
            y = Y/rMax; % in B(0,1)
            k1 = this.dilatation(rMax); % k1 must be approximated on {a < r < 1}
            
            % Radial quadrature (Bessel decomposition in 1D)
            rq = k1.radialQuadKernel(a,tol);
           
            
            % Circular quadrature 
            q2d = Quad2D(rq);
            
            % Local correction
            loc = localCorrections(x,y,a,k1,rq,tol,false);
            % We are ready to compute fast convolutions. 
            onlineEBD = @(v)(q2d.conv(x,y,v) + loc*v); 
        end
        function[MVx,MVy,rq,locx,locy] = offline_dEBD(this,X,Y,a,tol)
            
            % Rescaling
            rMax = rMaxCalc(X,Y); % diameter
            x = X/rMax; % in B(0,1)
            y = Y/rMax; % in B(0,1)
            k1 = this.dilatation(rMax); % k1 must be approximated on {a < r < 1}
            rq = k1.radialQuadKernel(a,tol*rMax,'grad',true);
            % the tolerance is set to tol*rMax because G'(rMax*u) = 1/rMax
            % k1'(u).
            q2d = Quad2D(rq);
            [q2dx,q2dy] = grad(q2d);
            loc = localCorrections(x,y,a,k1,rq,tol,true);
            locx = loc{1}; locy = loc{2};
            MVx = @(v)(-(q2dx.conv(x,y,v))/rMax + locx*v/rMax);
            MVy = @(v)(-(q2dy.conv(x,y,v))/rMax + locy*v/rMax);
        end
        
        % Display
        function[] = disp(this)
            fprintf('Kernel : function %s \n',func2str(this.func));            
        end
        
        function[out] = eval(this,x)
            fun = this.func;
            out = fun(x);
            out(abs(x) < 1e-12) = this.lim0;
        end
        function[out] = evalDer(this,x)
            fun = this.der;
            out = fun(x);
            out(abs(x) < 1e-12) = 0;
        end
        % Setters
        function[this] = setScalFunc(this,f)
            this.scalFunc = f;
        end
        function[this] = setNormFunc(this,g)
            this.normFunc = g;
        end
        function[this] = setStartFreq(this,k)
            this.startFreq = k;
        end
        % Addition of two kernels
        function[C] = plus(k1,k2)
            C = Kernel(@(x)(k1.func(x) + k2.func(x)),@(x)(k1.der(x) + k2.der(x)));
            sf1 = k1.scalFunc;
            sf2 = k2.scalFunc;
            C.scalFunc =  @(a,b,rho)(sf1(a,b,rho) + sf2(a,b,rho));
            C.customRadialQuad = @(a,tol,varargin)(k1.radialQuadKernel(a,tol/2,varargin{:}) ...
                + k2.radialQuadKernel(a,tol/2,varargin{:}));
            C.singular = or(k1.singular,k2.singular);
            C.lim0 = k1.lim0 + k2.lim0;
        end
        % Multiplication by a constant
        function[c] = mtimes(lambda,this)
            if and(isa(lambda,'double'),isscalar(lambda))
                assert(isa(this,'Kernel'));
                c = this;
                c.func = @(x)(lambda*this.func(x));
                c.der = @(x)(lambda*this.der(x));
                c.scalFunc = @(a,b,rho)(lambda*this.scalFunc(a,b,rho));
                c.normFunc = @(a,b)(abs(lambda)*this.normFunc(a,b));
                c.lim0 = lambda*this.lim0;
                c.customRadialQuad = @(a,tol,varargin)(lambda*this.radialQuadKernel(a,tol/lambda,varargin{:}));
            else
                assert(and(isa(this,'double'),isscalar(this)))
                c = times(this,lambda);
            end
        end
    end
    methods (Access = public)
        % Wrapper
        function[rq] = radialQuadKernel(this,a,tol,varargin)
            if ~ isempty(this.customRadialQuad)
                rq = this.customRadialQuad(a,tol,varargin{:});
            else
                rq = RadialQuadrature(a,this,tol,varargin{:});
            end
        end
        % Dilatation of the kernel : G(x) -> G(R*x).
        function[this] = dilatation(old,lambda)
            oldFunc = old.func;
            oldDer = old.der;
            ffunc = @(x)(oldFunc(lambda*x));
            dder = @(x)(lambda*oldDer(lambda*x));
            this = Kernel(ffunc,dder);
            this = this.setStartFreq(lambda*old.startFreq);
            this.lim0 = old.lim0;
        end
        
        
        
        
    end
end

