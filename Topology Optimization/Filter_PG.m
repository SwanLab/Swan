classdef Filter_PG < Filter
    properties
    end
    methods 
        function x_gp = getP0fromP1(obj,x)            
            x_gp = obj.P_operator*obj.Msmooth*x;
        end
        function x_reg = getP1fromP0(obj,x)
            x_reg = obj.P_operator'*obj.M0*x;
        end
    end
end