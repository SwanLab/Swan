classdef Filter_SLERP < Filter
    properties
    end
    methods 
        function x_gp = getP0fromP1(obj,x)
            M2=obj.faireF2(obj.coordinates',obj.connectivities',x);
            x_gp = obj.P_operator*M2;
        end
        function x_reg = getP1fromP0(obj,x)            
            x_reg = obj.P_operator'*obj.M0*x;
        end
        function x_reg= getP0fromP1_perimeter(obj,x,epsilon)
            obj.rhs = obj.faireF2(obj.coordinates',obj.connectivities',x);
            Rinv  = (epsilon^2*obj.Ksmooth + obj.Msmooth);
            x_reg = obj.solver.solve(Rinv,obj.rhs,obj.dof_per,obj.fixnodes_per);         
        end
    end
end