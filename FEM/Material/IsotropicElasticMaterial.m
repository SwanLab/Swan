classdef IsotropicElasticMaterial < ElasticMaterial

    properties (GetAccess = public, SetAccess = protected)
        kappa
        mu
        lambda
    end
    
    properties (Access = protected)
       E
       nu
       nstre        
    end
            
    methods (Access = public)        
        
        function obj = setProps(obj,props)
            obj.kappa = props.kappa;
            obj.mu    = props.mu;
            obj.computeLambda()
            obj.computeC();
        end
        
    end
    
    methods (Access = protected)
        
        function init(obj,cParams)
            obj.nelem = cParams.nelem;
            obj.nstre = cParams.nstre;
            obj.createCtensor();
        end
       
        function createCtensor(obj)
            obj.C = zeros(obj.nstre,obj.nstre,obj.nelem);
        end

    end
            
    methods (Access = protected, Abstract)
        computeC(obj)
        computeLambda(obj)
    end
    
end

