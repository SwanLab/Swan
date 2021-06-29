classdef Filter_PDE < Filter
    
    properties (Access = protected)
        Anodal2Gauss
        LHS
        epsilon
    end
    
    methods (Access = public)
        
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
           %!!!! EHHH 
           % x_reg = x;
           for igaus = 1:obj.quadrature.ngaus
                x0(:,igaus) = obj.Anodal2Gauss{igaus}*x_reg;                
           end
        end
        
        function x_reg = regularize(obj,F)
            RHS = obj.integrate_function_along_facets(F);
            x_reg = obj.solve_filter(RHS);
        end
        
        function obj = updateEpsilon(obj,epsilon)
            if obj.hasEpsilonChanged(epsilon)
                obj.epsilon = epsilon;
                obj.diffReacProb.setEpsilon(epsilon);
                obj.computeLHS();
            end
        end
        
    end
    
    methods (Access = protected)
        
        function computeLHS(obj)
            lhs = obj.diffReacProb.element.computeLHS();
            obj.LHS = decomposition(lhs);
        end                
        
    end
    
    methods (Access = private)
        
        function intX = integrate_P1_function_with_shape_function(obj,x)
            ndof = size(obj.Anodal2Gauss{1},2);
            intX = zeros(ndof,1);            
            for igaus = 1:obj.quadrature.ngaus
                dVG = obj.geometry.dvolu(:,igaus);
                xG = x(:,igaus);
                A = obj.Anodal2Gauss{igaus};
                intX = intX + A'*(xG.*dVG);
            end
        end
        
        function x_reg = solve_filter(obj,RHS)
            %x_reg = obj.LHS\(RHS);
            obj.diffReacProb.computeVariables(RHS);
            x_reg = obj.diffReacProb.variables.x;
        end

        
        function itHas = hasEpsilonChanged(obj,eps)
            if isempty(obj.epsilon)
                obj.epsilon = 0;
            end
            var = abs(eps - obj.epsilon)/eps;
            itHas = var > 1e-15;
        end        
        
    end
    
end