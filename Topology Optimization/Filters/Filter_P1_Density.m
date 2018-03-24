classdef Filter_P1_Density < Filter_P1
    properties
    end
    
    methods
        function obj = Filter_P1_Density(problemID,scale)
            obj@Filter_P1(problemID,scale);
        end
        
        function x_gp = getP0fromP1(obj,x)
            if isequal(x,obj.x)
                x_gp=obj.x_reg;
            else
                x_gp = obj.P_operator*obj.Msmooth*x;
                obj.x=x;
                obj.x_reg=x_gp;
            end
        end
        
        function x_reg = getP1fromP0(obj,x)
            gauss_sum=0;
            for igauss=1:length(x(1,:))
                gauss_sum=gauss_sum+obj.M0{igauss}*x(:,igauss);
            end
            x_reg = obj.P_operator'*gauss_sum;
        end
    end
end