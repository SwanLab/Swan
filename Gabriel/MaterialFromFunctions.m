classdef MaterialFromFunctions < handle

    properties (Access = private)
        density
        tensorFcn
        tensorDerivativeFcn
    end

    methods (Access = public)

        function obj = MaterialFromFunctions(cParams)
            obj.tensorFcn           = cParams.tensorFcn;
            obj.tensorDerivativeFcn = cParams.tensorDerivativeFcn;
        end

        function setDesignVariable(obj,x)
            obj.density = x{1};
        end

        function C = obtainTensor(obj)
            C = obj.tensorFcn(obj.density);
        end

        function dC = obtainTensorDerivative(obj)
            dC{1} = obj.tensorDerivativeFcn(obj.density);
        end

    end
end