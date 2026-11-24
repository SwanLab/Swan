classdef QuarterTurnParallelogramTransformation < handle

    methods (Access = public)

        function T = compute(~,params)
            a = params(1);
            s = params(2);
            turn  = [1-a,-a;a,1-a]/sqrt((1-a)^2+a^2);
            boost = [1,s;s,1]/sqrt(1-s^2);
            T = turn*boost;
        end

        function Tinv = computeInverse(obj,params)
            Tinv = inv(obj.compute(params));
        end

    end

end