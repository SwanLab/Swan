classdef Filter_PDE < Filter
    properties
        dvolu
        rhs
        A_nodal_2_gauss
    end
    
    methods
        function obj = Filter_PDE(problemID,scale)
            obj@Filter(problemID,scale);
        end
        
        function preProcess(obj)
            preProcess@Filter(obj);
            obj.dvolu = sparse(1:obj.diffReacProb.geometry.interpolation.nelem,1:obj.diffReacProb.geometry.interpolation.nelem,...
                sum(obj.diffReacProb.geometry.dvolu,2));
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
            rhs = (obj.A_nodal_2_gauss'*obj.M0{1}*x);
        end
        
        function x_reg = solve_filter(obj,rhs_x)
            obj.diffReacProb.computeVariables(rhs_x);
            x_reg = obj.diffReacProb.variables.x;
        end
        
        function obj = updateEpsilon(obj,epsilon)
            obj.diffReacProb.setEpsilon(epsilon);
        end
    end
end