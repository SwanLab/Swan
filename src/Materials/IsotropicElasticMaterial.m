classdef IsotropicElasticMaterial < Material
    
    properties (GetAccess = public, SetAccess = private)
        young
        poisson
        bulk
        shear
    end

    methods (Access = public)

        function obj = IsotropicElasticMaterial(cParams)
            obj.init(cParams)
            if isfield(cParams,'young') && isfield(cParams,'poisson')
                obj.young   = cParams.young;
                obj.poisson = cParams.poisson;
                obj.bulk  = obj.computeKappaFromYoungAndPoisson(obj.young,obj.poisson,obj.mesh.ndim);
                obj.shear = obj.computeMuFromYoungAndPoisson(obj.young,obj.poisson);
            elseif isfield(cParams,'bulk') && isfield(cParams,'shear')
                obj.bulk  = cParams.bulk;
                obj.shear = cParams.shear;
            else
                error('Young/Poisson or Bulk/Shear allowed')
            end
        end

    end

    methods (Access = protected)

        function C = evaluateNew(obj,xV)
            [l,m] = computeLameParameters(obj);
            lambda = l.evaluate(xV);
            mu = m.evaluate(xV);

            N = obj.mesh.ndim;
            nGauss = size(mu,2);
            nElem  = m.mesh.nelem;
            lambda = reshape(lambda,[1 1 1 1 nGauss nElem]);
            mu     = reshape(mu,[1 1 1 1 nGauss nElem]);
            I      = repmat(eye4D(N),[1 1 1 1 nGauss nElem]);
            IxI    = permute(repmat(kronEye(N),[1 1 1 1 nGauss nElem]),[4 3 2 1 5 6]);
            C = 2*mu.*I + lambda.*IxI;
        end

    end

    methods (Access = private)

        function [lambda,mu] = computeLameParameters(obj)
            lambda = obj.computeLambdaFromShearAndBulk(obj.shear,obj.bulk,obj.mesh.ndim);
            mu     = obj.shear;
        end

    end

    %Lame Parameters conversions
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

        function k = computeKappaFromShearAndLambda(m,l,N)
            k = (2/N)*m + l;
        end

        function mu = computeMuFromKappaAndLambda(k,l,N)
            mu = (N/2)*(k-l);
        end
        
    end
    
end

