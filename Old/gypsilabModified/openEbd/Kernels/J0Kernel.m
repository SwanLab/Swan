classdef J0Kernel < Kernel
    % Object of type kernel but optimized for J0 so that the computation of the 
    % radial quadrature goes faster
    % Radial quadrature is only made of 1 component !
    
    properties
        R=1, C=1;
        % Such that G(r) = C*besselj(0,R*r)
    end
    
    methods
        function[this] = J0Kernel(RR,CC)
            if nargin == 0
                RR = 1;
            end
            if nargin <= 1
                CC = 1;
            end
            this@Kernel(@(x)(CC*besselj(0,RR*x)),@(x)(-RR*CC*besselj(1,RR*x)))
            this.C = CC;
            this.R = RR;
        end
        function[out] = dilatation(this,lambda)
            out = J0Kernel(lambda*this.R,this.C);    
            % No changes in gamma_est, scalFunc nor normFunc
        end
        function[out] = mtimes(this,mu)
            if isa(this,'Kernel')
                out = J0Kernel(this.R,mu*this.C);
            else
                out = mtimes(mu,this);
            end
        end
        function[out] = radialQuadKernel(this,a,tol,varargin)
            Cmem = this.C;
            this = (1/Cmem)*this;
            rq = RadialQuadrature(a,this,tol,varargin{:});
            out = Cmem*rq;
        end
    end
    
    
end

