classdef Filter_SLERP_PDE < Filter_PDE
    properties
    end
    methods 
        function x_gp = getP0fromP1(obj,x)
            obj.rhs = obj.faireF2(obj.coordinates',obj.connectivities',x);
            Rinv  = (obj.epsilon^2*obj.Ksmooth + obj.Msmooth);
            dof.ndof=obj.dof_per.ndof;
            dof.vL=1:dof.ndof;
            dof.vR=[];
            fixnodes=[];
            x_reg = obj.solver.solve(Rinv,obj.rhs,dof,fixnodes);   
            x_gp=obj.A_nodal_2_gauss*x_reg;
        end
        function x_reg = getP1fromP0(obj,x)
            rhs = (obj.A_nodal_2_gauss'*obj.M0*x);
            Rinv = (obj.epsilon^2*obj.Ksmooth + obj.Msmooth);   
            dof.ndof=obj.dof_per.ndof;
            dof.vL=1:dof.ndof;
            dof.vR=[];
            fixnodes=[];
            x_reg = obj.solver.solve(Rinv,rhs,dof,fixnodes);    
        end
        function x_reg=getP0fromP1_per(obj,x,epsilon)
            obj.rhs = obj.faireF2(obj.coordinates',obj.connectivities',x);
            Rinv  = (epsilon^2*obj.Ksmooth + obj.Msmooth);
            x_reg = obj.solver.solve(Rinv,obj.rhs,obj.dof_per,obj.fixnodes_per);   
        end

    end
end