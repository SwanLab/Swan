classdef Y0Kernel < Kernel
    % Object of type kernel but optimized for Y0 so that the computation of the
    % radial quadrature goes faster
    
    properties (Access = public)
        C = 1;
        R; %such that G(r) = C*bessely(0,R*r)
    end
    
    methods
        function[kernel] = Y0Kernel(RR,CC)
            if nargin == 0
                RR = 1;
            end
            if nargin <= 1
                CC = 1;
            end
            func = @(x)(CC*bessely(0,RR*x));
            der = @(x)(-RR*CC*bessely(1,RR*x));
            kernel@Kernel(func,der);
            kernel.R = RR;
            kernel.C = CC;
            kernel = kernel.setScalFunc(@(a,b,rho)(CC*kernel.scalFuncHelmholtz(a,b,rho)));
            kernel = kernel.setNormFunc(@(a,b)(CC*kernel.normFuncHelmholtz(a,b)));
            kernel = kernel.setStartFreq(RR);
            kernel.singular = false;
            gamm = 0.577215664901532;
            kernel.lim0 = 2/pi*(gamm + log(RR/2)); 
        end
        
    end
    methods (Access = public)
        
        function[res] = scalFuncHelmholtz(this,a,b,rho)
            % Computes (a,b,\rho,R) -> \int_a^b x(Y_0(Rx))'J_1(\rho x) dx
            % ( = \int_a^b - Rx Y_1(Rx) J_1(\rho x) dx )
            RR = this.R;
            if ismember(RR,rho)
                error('bad choice of R ! It is not a root of Y0 or rho is not a root of J0')
            else
                res = helmholtzSP([a,b],rho,RR);
            end
        end
        function[res] = normFuncHelmholtz(this,a,b)
            % Computes (a,b) -> sqrt(2*pi*\int_a^b x [(Y0(Rx))']^2 dx))
            % (= sqrt(2*pi\int_a^b R^2 x Y_1(Rx)^2 dx)
            RR = this.R;
            T1 = b*RR*bessely(0,b*RR)^2;
            T2 = bessely(1,b*RR)*bessely(0,b*RR);
            T3 = b*RR*bessely(1,b*RR)^2;
            
            T4 = a*RR*bessely(0,a*RR)^2;
            T5 = bessely(1,a*RR)*bessely(0,a*RR);
            T6 = a*RR*bessely(1,a*RR)^2;
            
            res = sqrt(2*pi*RR/2*( b*(T1 - 2*T2 + T3) - a*(T4 - 2*T5 + T6)));
            
            % Explicit form
            %https://www.wolframalpha.com/input/?i=integral+from+a+to+b+of+R%5E2+x+*+(bessely(1,R+x))%5E2
        end
        function[out] = dilatation(this,lambda)
            lim0keep = this.lim0;
            out = Y0Kernel(lambda*this.R,this.C);
            out.lim0 = lim0keep;
        end
        function[out] = mtimes(this,mu)
            if isa(this,'Kernel')
                assert(and(isa(mu,'double'),isscalar(mu)));
                out = Y0Kernel(this.R,mu*this.C);
            else
                out = mtimes(mu,this);
            end
        end
        function[rq] = radialQuadKernel(this,a,tol,varargin)
            Cmem = this.C;
            this = (1/Cmem)*this;
            % 1 Detect if R is a root of Y0, and otherwise, find next root
            if abs(bessely(0,this.R))>1e-12
                if this.R < 0.89
                    c = -this.der(1)/this.func(1);
                    assert(c > 0);
                    rq = RadialQuadrature(a,this,tol,'robinCond',c,varargin{:});
                else
                    
                    Rtemp = nextY0root(this.R);
                    Ktemp = Y0Kernel(Rtemp);
                    % Dilatation factor
                    lambda = this.R/Rtemp;
                    atemp = a*lambda;
                    rqTemp = RadialQuadrature(atemp,Ktemp,tol/abs(Cmem),'startFreq',Rtemp,varargin{:});
                    % Use the previous radialQuad to form the final one :
                    rq = rqTemp.dilatation(lambda);
                end
            else
                rq = RadialQuadrature(a,this,tol/abs(Cmem),varargin{:},'startFreq',this.R);
            end
            rq = Cmem*rq;
        end
        
        
    end
end
