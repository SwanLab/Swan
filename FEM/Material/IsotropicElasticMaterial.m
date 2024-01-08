classdef IsotropicElasticMaterial < Material
    
    properties (SetAccess = private, GetAccess = public)
        young
        poisson
        bulk
        shear        
        lambda
    end

    properties (Access = protected)
        ndim                
        nstre
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
            obj.computeNstre();
            obj.computeOtherLameCoefficients();
        end

    end

    methods (Access = private)

        function computeNstre(obj)
            nDim      = obj.ndim;
            obj.nstre = 3*(nDim-1);
        end      

        function computeOtherLameCoefficients(obj)
            if isempty(obj.shear)
                obj.computeShear();
            end
            if isempty(obj.bulk)
                obj.computeBulk();
            end
            if isempty(obj.poisson)
                obj.computePoisson();
            end
            if isempty(obj.young)
                obj.computeYoung();
            end
            obj.computeLambda();
        end

        function computeShear(obj)
            E = obj.young;
            nu = obj.poisson;
            mu = obj.computeMuFromYoungAndPoisson(E,nu);
            obj.shear = mu;
        end        
        
        function computeBulk(obj)
            E  = obj.young;
            nu = obj.poisson;
            k  = obj.computeMuFromYoungAndPoisson(E,nu);
            obj.bulk = k;           
        end

        function computeYoung(obj)
            m = obj.shear;
            k = obj.bulk;
            E  = obj.computeYoungFromShearAndBulk(m,k);
            obj.young = E;           
        end

        function computePoisson(obj)
            m = obj.shear;
            k = obj.bulk;
            nu  = obj.computePoissonFromFromShearAndBulk(m,k);
            obj.poisson = nu;           
        end        

        function computeLambda(obj)
            m = obj.shear;
            k = obj.bulk;
            l  = obj.computeLambdaFromShearAndBulk(m,k);
            obj.lambda = l;           
        end        
       
        function mu = computeMuFromYoungAndPoisson(obj,E,nu)
            mu = E./(2*(1+nu));
        end

        function k = computeKappaFromYoungAndPoisson(obj,E,nu)
            N = obj.ndim;
            k = E./(N*(1-(N-1)*nu));
        end   

        function lambda = computeLambdaFromShearAndBulk(obj,m,k)
            N = obj.ndim;
            lambda = k - 2/N*m;
        end            

        function E = computeYoungFromShearAndBulk(obj,m,k)
            N = obj.ndim;            
            E = ((N*N*k).*(2*m))./(2*m + N*(N-1)*k);
        end
        
        function nu = computePoissonFromFromShearAndBulk(obj,m,k)
            N = obj.ndim;            
            nu = ((N*k)-(2*m))./(2*m + N*(N-1)*k);
        end     
        
    end
    
end

