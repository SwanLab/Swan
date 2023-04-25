classdef Element < handle 
    properties (Access = public)
        stifnessMatrix
    end
    properties(Access = protected)
        poissonCoefficient
        t
    end
    

    methods (Access = public)
        function obj = Element(cParams)
            obj.inputData(cParams);
        end
    end
    methods (Access = protected)
        function inputData(obj,cParams)
            obj.poissonCoefficient = cParams.poissonCoefficient; 
            obj.t = cParams.t; 
        end
    end
end