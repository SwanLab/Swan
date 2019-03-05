classdef Filter_PDE < Filter
    properties (Access = private)
        dvolu
        A_nodal_2_gauss
    end
    
    methods (Access = public)
        function preProcess(obj)
            preProcess@Filter(obj);
            obj.dvolu = sparse(1:obj.diffReacProb.geometry.interpolation.nelem,1:obj.diffReacProb.geometry.interpolation.nelem,...
                sum(obj.diffReacProb.geometry.dvolu,2));
            obj.A_nodal_2_gauss = obj.computeA;
        end
        
        function x_reg = getP1fromP1(obj,x)
            RHS = obj.integrate_L2_function_with_shape_function(x);
            x_reg = obj.solve_filter(RHS);
        end
        
        function x_reg = getP1fromP0(obj,x0)
            RHS = obj.integrate_P1_function_with_shape_function(x0);
            x_reg = obj.solve_filter(RHS);
        end
        
        function x0 = getP0fromP1(obj,x)
            x_reg =  obj.getP1fromP1(x);
            x0 = obj.A_nodal_2_gauss*x_reg;
        end
        
        function x_reg = regularize(obj,x,F)
            RHS = obj.integrate_function_along_facets(x,F);
            x_reg = obj.solve_filter(RHS);
        end
        
        function obj = updateEpsilon(obj,epsilon)
            obj.diffReacProb.setEpsilon(epsilon);
        end
    end
    
    methods (Access = private)
        function gauss_sum = integrate_P1_function_with_shape_function(obj,x)
            gauss_sum = 0;
            for igauss = 1:size(obj.M0,2)
                gauss_sum = gauss_sum+obj.A_nodal_2_gauss'*obj.M0{igauss}*x(:,igauss);
            end
        end
        
        function x_reg = solve_filter(obj,RHS)
            obj.diffReacProb.computeVariables(RHS);
            x_reg = obj.diffReacProb.variables.x;
        end
    end
end