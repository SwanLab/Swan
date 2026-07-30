classdef MultiMaterialInterpolator < handle

    methods (Static, Access = public)
        function [mu, kappa] = computeMuKappa(muRef,kRef,interp,N)
            nMat  = length(muRef);
            mu    = @(x) interp{nMat-2,nMat-1}.mu(x{end});
            kappa = @(x) interp{nMat-2,nMat-1}.k(x{end});
            for i = nMat-3:-1:1
                mu    = @(x) SimpAllInterpolator.computeMu(muRef{i},mu(x),kRef{i},kappa(x),x{i+1},N);
                kappa = @(x) SimpAllInterpolator.computeKappa(muRef{i},mu(x),kRef{i},kappa(x),x{i+1},N);
            end
            mu     = @(x) SimpAllInterpolator.computeMu(muRef{end},mu(x),kRef{end},kappa(x),x{1},N);
            kappa  = @(x) SimpAllInterpolator.computeKappa(muRef{end},mu(x),kRef{end},kappa(x),x{1},N);
        end

        % Improve with MuKappaDerivatives...

    end
end