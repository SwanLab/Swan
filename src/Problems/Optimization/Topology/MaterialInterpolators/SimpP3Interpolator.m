classdef SimpP3Interpolator < handle

    methods (Static, Access = public)

        function muI = computeMu(muA,muB,rho)
            muI = SimpP3Interpolator.interpolate(rho,muA,muB);
        end

        function kI = computeKappa(kA,kB,rho)
            kI = SimpP3Interpolator.interpolate(rho,kA,kB);
        end

        function dmuI = computeMuDerivative(muA,muB,rho)
            dmuI = SimpP3Interpolator.derive(rho,muA,muB);
        end

        function dkI = computeKappaDerivative(kA,kB,rho)
            dkI = SimpP3Interpolator.derive(rho,kA,kB);
        end

    end

    methods (Static, Access = private)
        function f = interpolate(rho,f0,f1)
            p = 3;
            [drho0,drho1] = SimpP3Interpolator.computeDensities(rho);
            f = (drho0^p)*f0 + (drho1^p)*f1;
        end

        function f = derive(rho,f0,f1)
            p = 3;
            [drho0,drho1] = SimpP3Interpolator.computeDensities(rho);
            f = -p*(drho0^(p-1))*f0 + p*(drho1^(p-1))*f1;
        end

        function [drho0,drho1] = computeDensities(rho)
            drho0 = 1 - rho;
            drho1 = rho;
        end
    end
end