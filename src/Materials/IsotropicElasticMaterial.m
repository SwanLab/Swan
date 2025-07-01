classdef IsotropicElasticMaterial < Material
    
    properties (SetAccess = private, GetAccess = private)
        young
        poisson
        bulk
        shear
    end

    properties (Access = protected)
        ndim
    end
    
    methods (Access = protected)

        function init(obj,cParams)
            obj.ndim    = cParams.ndim;
            if isfield(cParams,'young')
                obj.young = cParams.young;
            end
            if isfield(cParams,'poisson')
                obj.poisson = cParams.poisson;
            end            
            if isfield(cParams,'bulk')
                obj.bulk  = cParams.bulk;
            end
            if isfield(cParams,'shear')
                obj.shear = cParams.shear;
            end
        end

        function [mu,k] = computeShearAndBulk(obj,xV)
            if isempty(obj.shear) && isempty(obj.bulk)
                E  = obj.young.evaluate(xV);
                nu = obj.poisson.evaluate(xV);
                mu = obj.computeMuFromYoungAndPoisson(E,nu);
                k  = obj.computeKappaFromYoungAndPoisson(E,nu,obj.ndim);
            else
                mu = obj.shear.evaluate(xV);
                k  = obj.bulk.evaluate(xV); 
            end
        end


    end

    methods (Access = public, Static)
       
        function mu = computeMuFromYoungAndPoisson(E,nu)
            mu = E./(2*(1+nu));
        end

        function k = computeKappaFromYoungAndPoisson(E,nu,N)
            k = E./(N*(1-(N-1)*nu));
        end

        function E = computeYoungFromShearAndBulk(m,k,N)
            E = ((N*N*k).*(2*m))./(2*m + N*(N-1)*k);
        end
        
        function nu = computePoissonFromFromShearAndBulk(m,k,N)
            nu = ((N*k)-(2*m))./(2*m + N*(N-1)*k);
        end

        function lambda = computeLambdaFromShearAndBulk(m,k,N)
            lambda = k - 2/N*m;
        end
        
    end
    
end