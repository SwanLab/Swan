classdef Filter_P1_Density < Filter_P1
    
    methods (Access = public)
        
        function x_gp = getP0fromP1(obj,x)
            if isequal(x,obj.x)
                x_gp = obj.x_reg;
            else
                x_gp = obj.P_operator*obj.diffReacProb.element.M*x;
                obj.x = x;
                obj.x_reg = x_gp;
            end
        end
        
    end
    
end