classdef SimpAllInterpolator

    methods (Static, Access = public)

        function C = obtainTensor(matA,matB,rho,ndim)
            mesh = rho.mesh;
            op = @(xV) SimpAllInterpolator.evaluateTensor(matA,matB,rho,ndim,xV);
            C = DomainFunction.create(op,mesh,ndim^4);
        end

        function dC = obtainTensorDerivative(matA,matB,rho,ndim)
            mesh = rho.mesh;
            op = @(xV) SimpAllInterpolator.evaluateTensorDerivative(matA,matB,rho,ndim,xV);
            dC = DomainFunction.create(op,mesh,ndim^4);
        end

    end

    methods (Static, Access = private)

        function CxV = evaluateTensor(matA,matB,rho,N,xV)
            rhoV = rho.evaluate(xV);
            [mu,kappa] = SimpAllInterpolator.computeMuKappa(matA,matB,rhoV,N);
            CxV = SimpAllInterpolator.assembleTensor(mu,kappa,N);
        end

        function dCxV = evaluateTensorDerivative(matA,matB,rho,N,xV)
            rhoV = rho.evaluate(xV);
            [dmu,dkappa] = SimpAllInterpolator.computeMuKappaDerivative(matA,matB,rhoV,N);
            dCxV = SimpAllInterpolator.assembleTensor(dmu,dkappa,N);
        end

        
        function CxV = assembleTensor(mu,kappa,N)
            nGauss = size(mu,2);
            nElem  = size(mu,3);

            lambda = kappa - (2/N).*mu;

            I   = repmat(eye4D(N),  [1 1 1 1 nGauss nElem]);
            IxI = repmat(kronEye(N),[1 1 1 1 nGauss nElem]);

            muR     = reshape(mu,    [1 1 1 1 nGauss nElem]);
            lambdaR = reshape(lambda,[1 1 1 1 nGauss nElem]);

            CxV = 2*muR.*I + lambdaR.*IxI;
        end

        
        function [mu,kappa] = computeMuKappa(matA,matB,rho,N)
            m0 = matA.mu; k0 = matA.kappa;
            m1 = matB.mu; k1 = matB.kappa;

            eta0mu = SimpAllInterpolator.computeEtaMu(m0,k0,N);
            eta1mu = SimpAllInterpolator.computeEtaMu(m1,k1,N);
            cMu    = SimpAllInterpolator.computeCoeff(m0,m1,eta0mu,eta1mu);
            mu     = SimpAllInterpolator.computeRationalFunction(rho,cMu);

            eta0k  = SimpAllInterpolator.computeEtaKappa(m0,N);
            eta1k  = SimpAllInterpolator.computeEtaKappa(m1,N);
            cKappa = SimpAllInterpolator.computeCoeff(k0,k1,eta0k,eta1k);
            kappa  = SimpAllInterpolator.computeRationalFunction(rho,cKappa);
        end

        function [dmu,dkappa] = computeMuKappaDerivative(matA,matB,rho,N)
            m0 = matA.mu; k0 = matA.kappa;
            m1 = matB.mu; k1 = matB.kappa;

            eta0mu = SimpAllInterpolator.computeEtaMu(m0,k0,N);
            eta1mu = SimpAllInterpolator.computeEtaMu(m1,k1,N);
            cMu    = SimpAllInterpolator.computeCoeff(m0,m1,eta0mu,eta1mu);
            dmu    = SimpAllInterpolator.computeRationalDerivative(rho,cMu);

            eta0k  = SimpAllInterpolator.computeEtaKappa(m0,N);
            eta1k  = SimpAllInterpolator.computeEtaKappa(m1,N);
            cKappa = SimpAllInterpolator.computeCoeff(k0,k1,eta0k,eta1k);
            dkappa = SimpAllInterpolator.computeRationalDerivative(rho,cKappa);
        end

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