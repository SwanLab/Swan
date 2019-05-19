classdef Filter_P1_Density < Filter_P1
    
    properties (Access = private)
       Afilter 
    end
    
    methods (Access = public)
       
        function preProcess(obj)
           preProcess@Filter(obj)
           obj.Afilter = obj.P_operator*obj.diffReacProb.element.M;
        end
        
    end
    
    methods (Access = protected)
        
        function x0 = computeP0fromP1(obj,x)
            x0 = obj.Afilter*x;
            %x0 = obj.P_operator*obj.diffReacProb.element.M*x;
        end
        
    end
    
end