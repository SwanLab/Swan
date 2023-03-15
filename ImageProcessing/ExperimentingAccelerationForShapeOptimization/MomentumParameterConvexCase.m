classdef MomentumParameterConvexCase < MomentumParameter

    methods (Access = public, Static)

        function obj = MomentumParameterConvexCase(cParams)

        end

    end

    methods (Access = public)

        function v = computeValue(obj,cParams)
            iter = cParams.iter;
            v = iter/(iter+2);
        end

    end

end