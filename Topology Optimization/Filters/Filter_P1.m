classdef Filter_P1 < Filter
    
    methods (Access = protected, Abstract)
        computeP0fromP1(obj)
    end
    
    methods (Access = public)
        
        function obj = Filter_P1(cParams)
            obj@ Filter(cParams);
        end
        
        function preProcess(obj)
            preProcess@Filter(obj)
        end
        
        function x_reg = getP1fromP0(obj,x0)
            RHS = obj.integrate_P1_function_with_shape_function(x0);
            x_reg = obj.P_operator'*RHS;
        end
        
        function x0 = getP0fromP1(obj,x)
            if obj.xHasChanged(x)
                x0 = obj.computeP0fromP1(x);
            else
                x0 = obj.x_reg;
            end
            obj.updateStoredValues(x,x0);
        end
        
    end
    
    methods (Access = private)
        
        function intX = integrate_P1_function_with_shape_function(obj,x)
            intX = zeros(obj.nelem,1);
            for igaus = 1:obj.quadrature.ngaus
                dvolu = obj.geometry.dvolu(:,igaus);
                intX = intX + dvolu.*x(:,igaus);
            end
        end
        
    end
    
end