classdef H0Kernel < Kernel
    properties
        k; 
        C; % Such that G(r) = C*besselh(0,1,k*r)
        % Hankel function of the first kind. 
    end
    methods
        function[this] = H0Kernel(kk,CC)
            if nargin == 0
                kk = 1;
            end
            if nargin <= 1
                CC = 1;                                
            end
            this.k = kk;
            this.C = CC;
            Jk = J0Kernel(kk);
            Yk = Y0Kernel(kk);
            modelKern = this.C*(Jk + 1i*Yk);
            this.func = modelKern.func;
            this.der = modelKern.der;
            this.scalFunc = modelKern.scalFunc;
            this.lim0 = this.C*(Jk.lim0 + 1i*Yk.lim0);
        end
        function[out] = dilatation(this,lambda)
            lim0keep = this.lim0;
            out = H0Kernel(this.k*lambda,this.C);
            out.lim0 = lim0keep;
        end
        function[out] = mtimes(this,mu)
            if isa(this,'Kernel')
                assert(and(isa(mu,'double'),isscalar(mu)));
                out = H0Kernel(this.k,this.C*mu);
            else
                out = mtimes(mu,this);
            end
        end
        function[rq] = radialQuadKernel(this,a,tol,varargin)
            modelKern = J0Kernel(this.k) + 1i*Y0Kernel(this.k);
            rqTemp = modelKern.radialQuadKernel(a,tol/abs(this.C),varargin{:});
            rq = this.C*rqTemp;
        end
    end
end
