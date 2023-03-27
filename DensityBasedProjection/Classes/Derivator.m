classdef Derivator < handle
    properties (Access = public)
        derivedCost
    end
    properties (Access = private)
        derivingVariable
        derivatedVariable
        cParams
    end

    methods (Access = public)
        function obj = Derivator(cParams)
            obj.inputData(cParams);
        end

        function compute(obj)
            obj.derivate();
        end
    end
    methods (Access = private)
        function inputData(obj,cParams)
            obj.derivingVariable =  cParams.derivingVariable;
            obj.derivatedVariable =cParams.derivatedVariable;
            obj.cParams =cParams.cParams;
        end
        function derivate(obj)
            if obj.derivatedVariable == "ProjectedField" && obj.derivingVariable == "DerivedFilter"
                s.cParams = obj.cParams;
                B = CostFieldDerivator(s);
                B.compute();
                obj.derivedCost = B.derivedCost;
            elseif obj.derivatedVariable == "Cost" && obj.derivingVariable == "Field"
            elseif obj.derivatedVariable == "Cost" && obj.derivingVariable == "ProjectedField"
            elseif obj.derivatedVariable == "FilteredField" && obj.derivingVariable == "FieldDerivator"
            end
        end
    end
end