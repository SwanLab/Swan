classdef LameParametersConverter < handle

    methods (Static, Access = public)

        function kappa = computeBulkFromYoungAndPoisson(E,nu,N)
            kappa = E./(N*(1-(N-1)*nu));
        end

        function mu = computeShearFromYoungAndPoisson(E,nu)
            mu = E./(2*(1+nu));
        end

        function lambda = computeLambdaFromBulkAndShear(kappa,mu,N)
            lambda = kappa - (2/N)*mu;
        end

        function lambda = computeLambdaFromYoungAndPoisson(E,nu,N)
            kappa  = LameParametersConverter.computeBulkFromYoungAndPoisson(E,nu,N);
            mu     = LameParametersConverter.computeShearFromYoungAndPoisson(E,nu);
            lambda = LameParametersConverter.computeLambdaFromBulkAndShear(kappa,mu,N);
        end


    end
end