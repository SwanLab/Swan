classdef ComposedFunction < L2Function
      
    properties (Access = private)
        f1
        f2
        operation
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = ComposedFunction(cParams)
            obj.init(cParams)
        end

        function fxV = evaluate(obj, xG)
            f1V = obj.f1.evaluate(xG);
            f2V = obj.f2.evaluate(xG);
            fxV = obj.operation(f1V,f2V);
         end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.f1        = cParams.f1;
            obj.f2        = cParams.f2;
            obj.operation = cParams.operation;
        end
        
    end
    
end