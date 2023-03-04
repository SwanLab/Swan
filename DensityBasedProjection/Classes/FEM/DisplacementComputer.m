classdef DisplacementComputer < handle
    properties (Access = public)
        displacement
    end
    properties (Access = private)
        force
        globalStifnessMatrix 
        freeDegress
        elementNumberX
        elementNumberY
        
    end
    methods (Access = public)
        function obj = DisplacementComputer(cParams)
            obj.inputData(cParams);
        end
        function compute(obj)
            obj.computeDisplacement();
        end
    end
    methods (Access = private)
        function inputData(obj,cParams)
            obj.force = cParams.force;
            obj.globalStifnessMatrix =cParams.globalStifnessMatrix;
            obj.elementNumberX =cParams.elementNumberX;
            obj.elementNumberY =cParams.elementNumberY;
            obj.freeDegress =cParams.freeDegress;
        end
        function computeDisplacement(obj)
            obj.displacement = zeros(2*(obj.elementNumberY+1)*(obj.elementNumberX+1),1);
            obj.displacement(obj.freeDegress) = obj.globalStifnessMatrix(obj.freeDegress,obj.freeDegress)\obj.force(obj.freeDegress);
        end
    end
end