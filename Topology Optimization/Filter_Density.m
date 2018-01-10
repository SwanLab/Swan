classdef Filter_Density < Filter
    properties
    end
    methods 
        function x_gp = getP0fromP1(obj,x)     
            if norm(x) == norm(obj.x)
                x_gp=obj.x_reg;
            else
                x_gp = obj.P_operator*obj.Msmooth*x;
                obj.x=x;
                obj.x_reg=x_gp;
            end
        end
        function x_reg = getP1fromP0(obj,x)
            x_reg = obj.P_operator'*obj.M0*x;
        end
    end
end