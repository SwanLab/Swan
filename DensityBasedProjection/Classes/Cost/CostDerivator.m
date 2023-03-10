classdef CostDerivator < handle 
    properties (Access = public)
        derivedCost
    end
    properties (Access = private)
        derivedProyectedCost
        derivatedProjectedField
        derivedFilteredField
    end
    methods (Access = public)
        function obj = CostDerivator(cParams)
            obj.inputData(cParams);
        end

        function compute(obj)
            obj.derive();
        end
    end
    methods (Access = private)
        function inputData(obj,cParams)
            obj.derivedProyectedCost = cParams.derivedProyectedCost;
            obj.derivatedProjectedField = cParams.derivatedProjectedField;
            obj.derivedFilteredField = cParams.derivedFilteredField;
            obj.derivedCost = zeros(size(obj.derivedProyectedCost));
        end
        function derive(obj)
            obj.derivedCost(:) = obj.derivedFilteredField*(obj.derivedProyectedCost(:).*obj.derivatedProjectedField(:));
        end 
    end 
end