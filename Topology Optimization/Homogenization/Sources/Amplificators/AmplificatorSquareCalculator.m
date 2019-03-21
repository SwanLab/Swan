classdef AmplificatorSquareCalculator < handle
    
    properties (SetAccess = private, GetAccess = public)
        Phomog
        PhomogInv        
    end
    
    properties (Access = private)
       generalAmplificatorCalculator 
    end
    
    methods (Access = public)
        
        function obj = AmplificatorSquareCalculator(d)
            obj.init(d);
        end
        
        function compute(obj)
            obj.generalAmplificatorCalculator.compute();
            obj.transformAmplificatorInMatrixFormat();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,d)
            d.pNorm = 2;
            a = AmplificatorComponentsCalculator(d);
            obj.generalAmplificatorCalculator = a;            
        end
        
        function transformAmplificatorInMatrixFormat(obj)            
            Pv = obj.generalAmplificatorCalculator.Phomog;
            alpha = obj.generalAmplificatorCalculator.monom;
            for ialpha = 1:size(alpha,1)
                  component = alpha(ialpha,:)==1;
                  [i,j] = obj.indexVoigt2tensor(component);
                  Pm(i,j) = Pv(ialpha); 
                  Pm(j,i) = Pv(ialpha);                   
            end
            obj.Phomog = Pm;
            obj.PhomogInv = inv(Pm);
        end
        
    end
    
    methods (Access = private, Static)
        
        function [i,j] = indexVoigt2tensor(k)
            T = [1 1;
                2 2;
                3 3;
                2 3;
                1 3;
                1 2];
            i = T(k,1);
            j = T(k,2);
        end       
        
    end
end