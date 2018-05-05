classdef Filter_PDE < Filter
    properties
        solver
        rhs
        epsilon
        A_nodal_2_gauss
        element
    end
    methods
        function preProcess(obj,params)
            preProcess@Filter(obj,params);
            obj.P_operator=obj.computePoperator(obj.Msmooth);
            obj.element = params.element;
            obj.element.dof = params.dof;
            %obj.dof = physicalProblem.dof;
            obj.solver = Solver.create();
            obj.A_nodal_2_gauss = obj.computeA;
        end

        
        function x_reg = getP1fromP1(obj,x)
            rhs_x = obj.integrate_L2_function_with_shape_function(x);
            x_reg = obj.solve_filter(rhs_x);
        end
        
        function x_reg = getP1fromP0(obj,x)
            rhs_x = obj.integrate_P1_function_with_shape_function(x);
            x_reg = obj.solve_filter(rhs_x);
        end
        
        function x_gp = getP0fromP1(obj,x)
            x_reg =  obj.getP1fromP1(x);
            x_gp = obj.A_nodal_2_gauss*x_reg;
        end
        
        function rhs = integrate_P1_function_with_shape_function(obj,x)
            gauss_sum=0;
            for igauss=1:size(obj.M0,2)
                if size(x,2)==1
                    gauss_sum=gauss_sum+obj.A_nodal_2_gauss'*obj.M0{igauss}*x;
                else
                    gauss_sum=gauss_sum+obj.A_nodal_2_gauss'*obj.M0{igauss}*x(:,igauss);
                end
            end       
            rhs = gauss_sum;
        end
        
        % !! Can be done as a DiffReact_Problem !!
        function x_reg = solve_filter(obj,rhs_x)
            Rinv = (obj.epsilon^2*obj.Ksmooth + obj.Msmooth);
            Rinv_red = obj.element.full_matrix_2_reduced_matrix(Rinv);
            rhs_red  = obj.element.full_vector_2_reduced_vector(rhs_x);
            x_reg = obj.solver.solve(Rinv_red,rhs_red);
            x_reg = obj.element.reduced_vector_2_full_vector(x_reg);
        end        
    end
end