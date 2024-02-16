classdef MomentumParameterConstant < MomentumParameter

    properties (Access = private)
        value
    end

    methods (Access = public, Static)

        function obj = MomentumParameterConstant(cParams)
            obj.value = cParams.value;
        end

    end

    methods (Access = public)

        function v = computeValue(obj,cParams)
            v = obj.value;
        end
    end

end