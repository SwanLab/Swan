classdef IsotropicElasticMaterial < Material
    
    properties (Access = private)
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
            if isfield(cParams,'young') && isfield(cParams,'poisson')
                obj.young   = cParams.young;
                obj.poisson = cParams.poisson;
                obj.bulk  = obj.computeKappaFromYoungAndPoisson(obj.young,obj.poisson,obj.ndim);
                obj.shear = obj.computeMuFromYoungAndPoisson(obj.young,obj.poisson);
            elseif isfield(cParams,'bulk') && isfield(cParams,'shear')
                obj.bulk  = cParams.bulk;
                obj.shear = cParams.shear;
            else
                error('Young/Poisson or Bulk/Shear allowed')
            end
        end

        function [lambda,mu] = computeLameParameters(obj)
            lambda = obj.computeLambdaFromShearAndBulk(obj.shear,obj.bulk,obj.ndim);
            mu     = obj.shear;
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