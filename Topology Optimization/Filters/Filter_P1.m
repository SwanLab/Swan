classdef Filter_P1 < Filter
    properties
    end
    methods
        function preProcess(obj,physicalProblem)
            preProcess@Filter(obj,physicalProblem)
            obj.P_operator=obj.computePoperator(obj.Msmooth);
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