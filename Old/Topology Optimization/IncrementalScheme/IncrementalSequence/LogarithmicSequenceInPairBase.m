classdef LogarithmicSequenceInPairBase < IncrementalSequence
    
    properties (Access = private)
       coef 
    end
    
    methods (Access = public)
        
        function obj = LogarithmicSequenceInPairBase(a,nSteps,initialValue,finalValue)
            obj.coef = a;
            obj.nSteps = nSteps;
            obj.initialValue = initialValue;
            obj.finalValue = finalValue;            
            obj.generateAlphaSequence();
        end
        
        function update(obj,i)
            obj.value = (1-obj.alpha(i))*obj.initialValue + obj.alpha(i)*obj.finalValue;
            obj.value = round(obj.value/2)*2;
        end        
        
    end
    
    
    methods (Access = protected)
        
        function generateAlphaSequence(obj)
            %obj.alpha = logspace(obj.x0,obj.x1,obj.nSteps);
            %              exp = linspace(obj.x0,obj.x1,obj.nSteps);
            %              obj.alpha = 10.^exp;
            %              alpha0 = min(abs(obj.alpha));
            %              alpha1 = max(abs(obj.alpha));
            %              obj.alpha = (obj.alpha - alpha0)/(alpha1-alpha0);
            x = linspace(0,1,obj.nSteps);
            b = 1;
            a = b*(obj.coef+1)/(1-obj.coef);
            obj.alpha = (a*x./((a-1)*x +1)).^b;
        end
        
    end
    
end