classdef IsotropicElasticMaterial < Material
   
    properties (SetAccess = private, GetAccess = public)
        young
        poisson
        bulk
        shear 
    end

    properties (Access = protected)
        ndim
    end

    methods (Access = public)

        function mu = createShear(obj)
            s.operation = @(xV) obj.computeShear(xV);
            s.ndimf = 1;
            mu = DomainFunction(s);
        end

        function k = createBulk(obj)
            s.operation = @(xV) obj.computeBulk(xV);
            s.ndimf = 1;
            k = DomainFunction(s);
        end    

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

           

        function [muV,kV] = computeShearAndBulk(obj,xV)
            if isempty(obj.shear) && isempty(obj.bulk)
                mu = obj.createShear();
                k  = obj.createBulk();
            else
                mu = obj.shear;
                k  = obj.bulk;
            end
            muV = mu.evaluate(xV);
            kV  = k.evaluate(xV);
        end


    end

    methods (Access = private)

        function mu = computeShear(obj,xV)
            E  = obj.young.evaluate(xV);
            nu = obj.poisson.evaluate(xV);
            mu = obj.computeMuFromYoungAndPoisson(E,nu);
        end

        function k = computeBulk(obj,xV)
            E  = obj.young.evaluate(xV);
            nu = obj.poisson.evaluate(xV);
            N  = obj.ndim;
            k  = obj.computeKappaFromYoungAndPoisson(E,nu,N);
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

