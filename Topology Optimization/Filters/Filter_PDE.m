classdef Filter_PDE < Filter
    
    properties (Access = private)
        Anodal2Gauss
    end
    
    methods (Access = public)
        
        function obj = Filter_PDE(cParams)
            obj@Filter(cParams);
        end
        
        function preProcess(obj)
            preProcess@Filter(obj);
            obj.Anodal2Gauss = obj.computeA();
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
            x0 = obj.Anodal2Gauss*x_reg;
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
        
        function intX = integrate_P1_function_with_shape_function(obj,x)
            ndof = size(obj.Anodal2Gauss,2);
            intX = zeros(ndof,1);            
            for igaus = 1:obj.quadrature.ngaus
                dv = obj.geometry.dvolu(:,igaus);
                intX = intX + obj.Anodal2Gauss'*(dv.*x(:,igaus));
            end
        end
        
        function x_reg = solve_filter(obj,RHS)
            obj.diffReacProb.computeVariables(RHS);
            x_reg = obj.diffReacProb.variables.x;
        end
        
    end
    
end