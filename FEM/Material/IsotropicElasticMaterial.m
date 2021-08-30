classdef IsotropicElasticMaterial < ElasticMaterial

    properties (GetAccess = public, SetAccess = protected)
        nstre                       
    end
    
    properties (Access = protected)
        kappa
        mu
        lambda       
    end
            
    methods (Access = public)        
        
        function compute(obj,s)
            obj.kappa = s.kappa;
            obj.mu    = s.mu;
            obj.nElem = size(obj.mu,1);
            obj.nGaus = size(obj.mu,2);            
            obj.computeC();
        end
        
    end
    
    methods (Access = protected)
        
        function init(obj,cParams)
            obj.nstre = cParams.nstre;
        end
        
    end
            
    methods (Access = protected, Abstract)
        computeC(obj)
    end
    
    methods (Access = public, Static)
       
        function mu = computeMuFromYoungAndNu(E,nu)
            mu = E./(2*(1+nu));
        end        
    end
    
end

