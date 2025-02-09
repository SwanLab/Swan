classdef ThinPlate < Kernel
    % Object of type kernel but optimized for ThinPlate so that the computation of the
    % radial quadrature goes faster
    properties (Access = public)
        C = 1;  %such that G(r) = C*r^2*log(R*r)
        R = 1;
    end
    
    methods
        function[kernel] = ThinPlate(CC,RR)
            if nargin == 0
                CC = 1;
            end
            if nargin <= 1
                RR = 1;
            end
            kernel@Kernel(@(x)(CC*x.^2.*log(RR*x)),@(x)(2*CC*x.*log(RR*x) + CC*x))
            kernel.scalFunc = @(a,b,rho)(CC*thinPlateSP([a,b],rho,RR));
            kernel.normFunc = @(a,b)(abs(CC)*ThinPlate.norm_func(a,b,RR));
            kernel.gamma_est = @ThinPlate.gamma_est;
            kernel.C = CC;
            kernel.R = RR;
            kernel.singular = false;
        end
    end
    methods (Access = public)
        function[out] = dilatation(this,lambda)
            CC = this.C;
            RR = this.R;
            out = ThinPlate(lambda^2*CC,lambda*RR);
        end
        function[out] = mtimes(this,mu)
            if isa(this,'Kernel')
                assert(and(isa(mu,'double'),isscalar(mu)));
                out = ThinPlate(mu*this.C);
            else
                out = mtimes(mu,this);
            end
        end
        function[out] = radialQuadKernel(this,a,tol,varargin)
            Cmem = this.C;
            this = (1/Cmem)*this;
            rq = RadialQuadrature(a,this,tol/abs(Cmem),varargin{:});
            out = Cmem*rq;
        end
        function[onlineEBD,rq,loc] = offlineEBD(this,X,Y,a,tol)
            % Special way to handle the radial quadrature : we add a term
            % in r^2 to enforce multi-Dirichlet condition.
            rMax = rMaxCalc(X,Y);
            x = X/rMax;
            y = Y/rMax;
            k1 = ThinPlate(this.C*rMax^2,this.R*rMax) +...
                X2kKernel(1,-this.C*(1+log(rMax*this.R))*rMax^2); % adjust multi-Dirichlet condition.
            rq = RadialQuadrature(a,k1,tol);
            q2d = Quad2D(rq);
            loc = localCorrections(x,y,a,k1,rq,tol,false);
            onlineEBD = @(v)(q2d.conv(x,y,v) + loc*v + suppTerm(v));
            
            % We remove the added r^2 term by computing the convolution as
            % follows:
            function[res] = suppTerm(alpha)
                alphay = [alpha alpha].*y;
                x2 = x(:,1).^2 + x(:,2).^2;
                y2 = y(:,1).^2 + y(:,2).^2;
                alphay2 = alpha.*y2;
                res = this.C*rMax^2*(log(rMax*this.R)+1)*(...
                    x2*sum(alpha)...
                    - 2*x(:,1)*sum(alphay(:,1)) - 2*x(:,2)*sum(alphay(:,2))...
                    + sum(alphay2));
            end
        end
    end
    methods (Static, Access = protected)
        function[res] = norm_func(a,b,R)
            % Found with Maple
            if a==0
                res = sqrt(2*pi*(...
                    b ^ 4 * log(R * b) ^ 2 +...
                    b ^ 4 * log(R * b) / 0.2e1 + b ^ 4 / 0.8e1));
                
            else
                res = sqrt(2*pi*(...
                    -a ^ 4 * log(R * a) ^ 2 - a ^ 4 * log(R * a) / 0.2e1 ...
                    - a ^ 4 / 0.8e1 + b ^ 4 * log(R * b) ^ 2 +...
                    b ^ 4 * log(R * b) / 0.2e1 + b ^ 4 / 0.8e1));
            end
        end
        function[low,up] = gamma_est(~)
            % Helps the radial quadrature to guess the number of components
            up = 7;
            low= 0;
        end
        
    end
    
    
end

