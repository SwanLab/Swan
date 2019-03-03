classdef Filter_P1_Density < Filter_P1
    
    methods (Access = protected)
        
        function x0 = computeP0fromP1(obj,x)
            x0 = obj.P_operator*obj.diffReacProb.element.M*x;
        end
        
    end
    
end