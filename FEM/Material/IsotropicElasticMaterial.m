classdef IsotropicElasticMaterial < ElasticMaterial

    properties (GetAccess = public, SetAccess = protected)
       nstre        
    end
    
    properties (Access = protected)
        ndim
       
        kappa
        mu
        lambda       
    end
            
    methods (Access = public)        
        
        function compute(obj,s)
            obj.kappa = s.kappa;
            obj.mu    = s.mu;
            obj.computeLambda();
            obj.computeC(obj.mu,obj.lambda);
        end
        
    end
    
    methods (Access = protected)
        
        function init(obj,cParams)
            obj.nElem = cParams.nelem;
            obj.nstre = cParams.nstre;
        end
       
    end
            
    methods (Access = protected, Abstract)
        computeC(obj)
        computeLambda(obj)
    end
    
end

