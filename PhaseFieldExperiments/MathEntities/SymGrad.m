classdef SymGrad < handle

    properties (Access = private)
        F
    end
    
    methods (Access = public)
        
        function obj = SymGrad(F)
            obj.init(F)
            
        end
        
        function symGradF = evaluate(obj,xV)
            symGradF = 0.5*(Grad(obj.F).evaluate(xV) + pagetranspose(Grad(obj.F).evaluate(xV)));
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,F)
            obj.F = F;
        end
        
    end
    
end