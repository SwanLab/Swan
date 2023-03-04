classdef Cost < handle
    properties (Access = protected)
        force
        displacement
        globalStifnessMatrix
    end
    methods (Access = public)
        function obj = Cost(cParams)
            obj.inputData(cParams);
        end
    end
    methods (Access = private)
        function inputData(obj,cParams)
            obj.force = cParams.force;
            obj.displacement =cParams.displacement;
        end
    end
end