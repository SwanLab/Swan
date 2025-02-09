classdef Constraint2 < CC
    
    properties (Access = public)
        designVariable
        nElem
    end
    
    methods (Access = public)
        
        function obj = Constraint2(cParams)
            obj.init(cParams)
        end
        
        function computeFunctionAndGradient(obj,iter)
            obj.computeFunctions(iter);
            obj.computeGradients(iter);
        end
        
    end
    
end