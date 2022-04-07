classdef Cost < CC
    
    properties (Access = private)
        designVariable
        nElem
    end
    
    methods (Access = public)
        
        function obj = Cost(cParams)
            obj.init(cParams)
        end
        
        function computeFunctionAndGradient(obj,settings)
            obj.computeFunctions(settings);
            obj.computeGradients(settings);
        end

    end
    
end