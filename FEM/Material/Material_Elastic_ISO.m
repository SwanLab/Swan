classdef Material_Elastic_ISO < Material_Elastic

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
            obj.mu = props.mu;
            obj.lambda = obj.kappa-obj.mu;
            obj.computeC();
        end
        
    end
    
    methods (Access = protected)
       
        function createCtensor(obj)
            obj.C = zeros(obj.nstre,obj.nstre,obj.nelem);
        end

    end
            
    methods (Access = protected, Abstract)
        computeC(obj)
    end
    
end

