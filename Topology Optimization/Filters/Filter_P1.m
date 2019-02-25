classdef Filter_P1 < Filter
    
    methods (Access = public)
        
        function preProcess(obj)
            preProcess@Filter(obj)
            obj.P_operator = obj.computePoperator(obj.diffReacProb.element.M);
        end
        
        function x_reg = getP1fromP0(obj,x)
            RHS = obj.integrate_P1_function_with_shape_function(x);
            x_reg = obj.P_operator'*RHS;
        end
        
    end
    
    methods (Access = private)
        
        function gauss_sum = integrate_P1_function_with_shape_function(obj,x)
            gauss_sum=0;
            for igauss=1:size(obj.M0,2)
                gauss_sum=gauss_sum+obj.M0{igauss}*x(:,igauss);
            end
        end
        
    end
    
end