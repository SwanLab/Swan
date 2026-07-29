classdef SimpAllInterpolator < handle

    methods (Static, Access = public)

        function muI = computeMu(muA,muB,kA,kB,rho,N)
            eta0mu = SimpAllInterpolator.computeEtaMu(muA,kA,N);
            eta1mu = SimpAllInterpolator.computeEtaMu(muB,kB,N);
            cMu    = SimpAllInterpolator.computeCoeff(muA,muB,eta0mu,eta1mu);
            muI    = SimpAllInterpolator.computeRationalFunction(rho,cMu);
        end

        function kI = computeKappa(muA,muB,kA,kB,rho,N)
            eta0k  = SimpAllInterpolator.computeEtaKappa(muA,N);
            eta1k  = SimpAllInterpolator.computeEtaKappa(muB,N);
            cKappa = SimpAllInterpolator.computeCoeff(kA,kB,eta0k,eta1k);
            kI     = SimpAllInterpolator.computeRationalFunction(rho,cKappa);
        end

        function dmuI = computeMuDerivative(muA,muB,kA,kB,rho,N)
            eta0mu = SimpAllInterpolator.computeEtaMu(muA,kA,N);
            eta1mu = SimpAllInterpolator.computeEtaMu(muB,kB,N);
            cMu    = SimpAllInterpolator.computeCoeff(muA,muB,eta0mu,eta1mu);
            dmuI   = SimpAllInterpolator.computeRationalDerivative(rho,cMu);
        end

        function dkI = computeKappaDerivative(muA,muB,kA,kB,rho,N)
            eta0k  = SimpAllInterpolator.computeEtaKappa(muA,N);
            eta1k  = SimpAllInterpolator.computeEtaKappa(muB,N);
            cKappa = SimpAllInterpolator.computeCoeff(kA,kB,eta0k,eta1k);
            dkI    = SimpAllInterpolator.computeRationalDerivative(rho,cKappa);
        end

    end

    methods (Static, Access = private)

        function etaMu = computeEtaMu(mu,kappa,N)
            num = -mu.*(4.*mu - kappa.*N^2 - 2.*mu.*N^2 + 2.*mu.*N);
            den = 2.*N.*(kappa + 2*mu);
            etaMu = num./den;
        end

        function etaKappa = computeEtaKappa(mu,N)
            etaKappa = 2.*mu.*(N-1)./N;
        end

        function c = computeCoeff(f0,f1,eta0,eta1)
            c.n01 = -(f1 - f0).*(eta1 - eta0);
            c.n0  = f0.*(f1 + eta0);
            c.n1  = f1.*(f0 + eta1);
            c.d0  = (f1 + eta0);
            c.d1  = (f0 + eta1);
        end

        function f = computeRationalFunction(rho,c)
            num = c.n01.*(1-rho).*rho + c.n0.*(1-rho) + c.n1.*rho;
            den = c.d0.*(1-rho) + c.d1.*rho;
            f   = num./den;
        end

        function df = computeRationalDerivative(rho,c)
            num    = c.n01.*(1-rho).*rho + c.n0.*(1-rho) + c.n1.*rho;
            den    = c.d0.*(1-rho) + c.d1.*rho;
            derNum = (c.n01.*(1-2.*rho)-c.n0+c.n1).*den - num.*(-c.d0+c.d1);
            df     = derNum./(den.^2);
        end

    end
end