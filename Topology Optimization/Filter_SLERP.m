classdef Filter_SLERP < Filter
    properties
    end
    methods 
        function x_gp = getP0fromP1(obj,x)
            if norm(x) == norm(obj.x)
                x_gp=obj.x_reg;
            else
                M2=obj.faireF2(obj.coordinates',obj.connectivities',x);
                x_gp = obj.P_operator*M2;
                obj.x=x;
                obj.x_reg=x_gp;
            end
        end
        function x_reg = getP1fromP0(obj,x)
            
            x_reg = obj.P_operator'*obj.M0{1}*x;
        end
        function x_reg= getP0fromP1_perimeter(obj,x,epsilon)
            obj.rhs = obj.faireF2(obj.coordinates',obj.connectivities',x);
            Rinv  = (epsilon^2*obj.Ksmooth + obj.Msmooth);
            x_reg = obj.solver.solve(Rinv,obj.rhs,obj.dof_per,obj.fixnodes_per);         
        end
    end
end