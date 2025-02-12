classdef X2kKernel < Kernel
    % Object of type kernel but optimized for log so that the computation of the 
    % radial quadrature goes faster
    properties (Access = public)
        C = 1; 
        k = 1; % such that G(r) =C*r^(2*k) 
    end
    
    methods
      	function[kernel] = X2kKernel(kk,CC)
            if nargin == 0
                kk = 1;
            end
            if nargin <= 1
                CC = 1;
            end
            kernel@Kernel(@(x)(CC*x.^(2*kk)),@(x)(2*kk*CC*x.^(2*kk - 1)))
            kernel.scalFunc = @(a,b,rho)(CC*X2kSP([a,b],rho,kk));
            kernel.normFunc = @(a,b)(abs(CC)*sqrt(2*pi*4*kk^2/(4*kk - 1)*(b^(4*kk - 1) - a^(4*kk- 1))));
            kernel.gamma_est = @X2kKernel.gamma_est;
            kernel.k = kk;
            kernel.C = CC;
            kernel.singular = false;
        end
    end
    methods (Access = public)
        function[out] = dilatation(this,lambda)
            out = X2kKernel(this.k,lambda^(2*this.k));
            % No changes in gamma_est, scalFunc nor normFunc
        end
        function[out] = mtimes(this,mu)
            if isa(this,'Kernel')
                assert(and(isa(mu,'double'),isscalar(mu)));
                out = X2kKernel(this.k,mu*this.C);
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
        function[low,up] = gamma_est(~)
            % Helps the radial quadrature to guess the number of components
            up = 4;
            low= 0;
        end
    end
    
end


