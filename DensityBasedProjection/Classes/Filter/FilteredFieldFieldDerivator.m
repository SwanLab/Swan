classdef FilteredFieldFieldDerivator < handle 
    properties (Access = public)
        derivedFilteredField
    end
    properties (Access = private)
        H
        Hs        
    end
    methods (Access = public)
        function obj = FilteredFieldFieldDerivator(cParams)
            obj.inputData(cParams);
        end

        function compute(obj)
            obj.derive();
        end
    end
    methods (Access = private)
        function inputData(obj,cParams)
            obj.H = cParams.H;
            obj.Hs = cParams.Hs;
        end
        function derive(obj)
            obj.derivedFilteredField = (obj.H./(obj.Hs'));
        end 
    end 
end