classdef LineSearchInitiator < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = protected)
       designVariable
       objectiveFunction       
    end
    
    methods (Access = public, Static)
        
        function obj = create(cParams)
           f = LineSearchInitiatorFactory();
           obj = f.create(cParams);
        end
        
    end
    
    methods (Access = protected)
        
        function init(obj,cParams)
            obj.designVariable    = cParams.designVariable;
            obj.objectiveFunction = cParams.objectiveFunction;            
        end
        
    end
    
end