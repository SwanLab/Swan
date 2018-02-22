classdef Filter_PDE < Filter
    properties
        dof_per
        solver
        rhs
        epsilon
        A_nodal_2_gauss
        element
    end
    methods
        function preProcess(obj,physicalProblem)
            preProcess@Filter(obj,physicalProblem);
            obj.element = physicalProblem.element;
            obj.dof_per=physicalProblem.dof;
            %obj.dof = physicalProblem.dof;
            obj.solver = Solver.create();
            obj.epsilon=0.03;
            obj.A_nodal_2_gauss=obj.computeA(physicalProblem);
        end

        
        function x_reg=getP1fromP1(obj,x)
            rhs_x = obj.integrate_L2_function_with_shape_function(x);
            x_reg = obj.solve_filter(rhs_x);
        end
        
        function x_reg = getP1fromP0(obj,x)
            rhs_x = obj.integrate_P1_function_with_shape_function(x);
            x_reg = obj.solve_filter(rhs_x);
        end
        
        function x_gp = getP0fromP1(obj,x)
            x_reg= obj.getP1fromP1(x);
            x_gp = obj.A_nodal_2_gauss*x_reg;
        end
        
        function rhs = integrate_P1_function_with_shape_function(obj,x)
            rhs = (obj.A_nodal_2_gauss'*obj.M0{1}*x);
        end
        
        function x_reg = solve_filter(obj,rhs_x)
            Rinv = (obj.epsilon^2*obj.Ksmooth + obj.Msmooth);
            Rinv_red = obj.element.full_matrix_2_reduced_matrix(Rinv,obj.dof_per);
            rhs_red  = obj.element.full_vector_2_reduced_vector(rhs_x,obj.dof_per);
            x_reg = obj.solver.solve(Rinv_red,rhs_red);
            x_reg = obj.element.reduced_vector_2_full_vector(x_reg,obj.dof_per);
        end
        

        
    end
end