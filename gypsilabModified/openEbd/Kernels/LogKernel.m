classdef LogKernel < Kernel
    % Object of type kernel but optimized for log so that the computation of the 
    % radial quadrature goes faster
    properties (Access = public)
        C = 1;
        R = 1;  %such that G(r) =  C*log(R*r) 
    end
    
    methods
      	function[kernel] = LogKernel(RR,CC)
            if nargin == 0
                RR = 1;
            end
            if nargin <= 1
                CC = 1;
            end
            kernel@Kernel(@(x)(CC*log(RR*x)),@(x)(CC*1./x))
            kernel.scalFunc = @(a,b,rho)(CC*laplaceSP([a,b],rho));
            kernel.normFunc = @(a,b)(abs(CC)*sqrt(2*pi*log(b/a)));
            kernel.gamma_est = @LogKernel.gamma_est;
            kernel.R = RR;
            kernel.C = CC;
            kernel.singular = false;
        end
    end
    methods (Access = public)
        function[out] = dilatation(this,lambda)
            out = LogKernel(lambda*this.R,this.C);    
            % No changes in gamma_est, scalFunc nor normFunc
        end
        function[out] = mtimes(this,mu)
            if isa(this,'Kernel')
                assert(and(isa(mu,'double'),isscalar(mu)));
                out = LogKernel(this.R,mu*this.C);
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
    end
    methods (Access = protected)
         
    end
    methods (Static, Access = protected)
        function[low,up] = gamma_est(tol)
            % Helps the radial quadrature to guess the number of components
            up = 1/5.8*log(180/tol);
            low= max(up-1,0.5);
        end
    end
    
end

