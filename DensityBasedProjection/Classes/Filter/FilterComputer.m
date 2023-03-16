classdef FilterComputer < handle 
    properties (Access = public)
        filteredField
    end
    properties (Access = private)
        filterParameters
        field
    end
    methods (Access = public)
        function obj = FilterComputer(cParams)
            obj.inputData(cParams);
        end

        function compute(obj)
            obj.filter();
        end
    end
    methods (Access = private)
        function inputData(obj,cParams)
            obj.filterParameters =cParams.filterParameters;
            obj.field =cParams.field;
            obj.filteredField = zeros(size(obj.field));
        end
        function filter(obj)
            obj.filteredField(:) = (obj.filterParameters.H*obj.field(:))./obj.filterParameters.Hs;
        end 
    end
end