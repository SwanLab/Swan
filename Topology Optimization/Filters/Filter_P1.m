classdef Filter_P1 < Filter
    properties
        P_operator
    end
    
    methods
        function obj = Filter_P1(problemID,scale)
            obj@Filter(problemID,scale);
        end
        
        function preProcess(obj)
            preProcess@Filter(obj)
            obj.P_operator = obj.computePoperator(obj.diffReacProb.element.M);
        end
        function x_reg = getP1fromP0(obj,x)
            gauss_sum=0;
            for igauss=1:size(obj.M0,2)
                if size(x,2)==1
                    gauss_sum=gauss_sum+obj.M0{igauss}*x;
                else
                    gauss_sum=gauss_sum+obj.M0{igauss}*x(:,igauss);
                end
            end
            x_reg = obj.P_operator'*gauss_sum;
        end
        
    end
end