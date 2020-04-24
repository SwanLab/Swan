classdef Momentum < handle
    
    properties (Access = private)
        iterator
        designVariable
        momentumParameter
        tOld
        tOldOld
    end
    
    methods (Access = public)
        
        function obj = Momentum(cParams)
            obj.designVariable = cParams.designVariable;
            obj.iterator       = cParams.iterator;
        end
        
        function apply(obj)
            obj.computeMomentumParameter();
            obj.addIntertialTerm();
        end
        
    end
    
    methods (Access = private)
        
        function computeMomentumParameter(obj)
            i = obj.iterator.value;
            if i < 3
                obj.tOld = 1;
                obj.tOldOld = 1;
                beta = 1;
            else
                t = (1+sqrt(1+4*obj.tOld^2))/2;
                beta = (obj.tOldOld-1)/t;  % beta = i/(i+3);
                obj.tOldOld = obj.tOld;
                obj.tOld = t;
            end
            obj.momentumParameter = beta;
        end
        
        function addIntertialTerm(obj)
            beta    = obj.momentumParameter;
            xOldOld = obj.designVariable.valueOldOld;
            xOld    = obj.designVariable.valueOld;
            x = xOld + beta*(xOld - xOldOld);
            obj.designVariable.value = x;
        end
        
    end
    
    
    
end