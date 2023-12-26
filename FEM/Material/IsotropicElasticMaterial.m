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
            obj.young   = cParams.young;
            obj.poisson = cParams.poisson; 
            obj.computeNstre();
            obj.computeShear();
            obj.computeBulk();     
            obj.computeLambda();
        end

    end

    methods (Access = private)

        function computeNstre(obj)
            nDim      = obj.ndim;
            obj.nstre = 3*(nDim-1);
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
            E = ((N*N*k)*(2*m))./(2*m + N*(N-1)*k);
        end
        
        function nu = computeNuFromFromShearAndBulk(obj,m,k)
            N = obj.ndim;            
            nu = ((N*k)-(2*m))./(2*m + N*(N-1)*k);
        end     
        
    end
    
end

