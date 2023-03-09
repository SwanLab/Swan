classdef ForceComputer < handle
    properties (Access = public)
        force
    end
    properties (Access = private)
        neumanCondition
        output
        elementNumberX
        elementNumberY
    end
    methods (Access = public)
        function obj = ForceComputer(cParams)
            obj.inputData(cParams);
        end
        function compute(obj)
            obj.computeForce();
        end
    end
    methods (Access = private)
        function inputData(obj,cParams)
            obj.neumanCondition = cParams.neumanCondition;
            obj.output =cParams.output;
            obj.elementNumberX =cParams.elementNumberX;
            obj.elementNumberY =cParams.elementNumberY;
        end
        function computeForce(obj)
            obj.force  = sparse(obj.output,1,obj.neumanCondition,2*(obj.elementNumberY+1)*(obj.elementNumberX+1),1);
        end
    end
end