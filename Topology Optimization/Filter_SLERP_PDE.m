classdef Filter_SLERP_PDE < Filter_PDE
    properties
    end
    methods
        function x_gp = getP0fromP1(obj,x,epsilon)
            obj.rhs = obj.faireF2(obj.coordinates',obj.connectivities',x);
%             Rinv  = (epsilon^2*obj.Ksmooth + obj.Msmooth);
             x_reg = zeros(obj.dof_per.ndof,1);
%             x_reg = obj.solver.solve(x_reg,Rinv,obj.rhs,obj.dof_per);
%             x_gp=obj.A_nodal_2_gauss*x_reg;
        end
        function x_reg = getP1fromP0(obj,x)
            rhs = (obj.A_nodal_2_gauss'*obj.M0*x);
            Rinv = (obj.epsilon^2*obj.Ksmooth + obj.Msmooth);
            dof.ndof=obj.dof_per.ndof;
            dof.vL=1:dof.ndof;
            dof.vR=[];
            x_reg = zeros(obj.dof.ndof,1);
            x_reg = obj.solver.solve(x_reg,Rinv,rhs,dof);
        end
        function x_reg=getP0fromP1_per(obj,x,epsilon)
            obj.rhs = obj.faireF2(obj.coordinates',obj.connectivities',x);
            Rinv  = (epsilon^2*obj.Ksmooth + obj.Msmooth);
            x_reg = zeros(obj.dof_per.ndof,1);
            x_reg = obj.solver.solve(x_reg,Rinv,obj.rhs,obj.dof_per);
        end
        
    end
end