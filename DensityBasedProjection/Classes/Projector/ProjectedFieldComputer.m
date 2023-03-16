classdef ProjectedFieldComputer < handle
    properties (Access = protected)
        beta 
        eta
        filteredField      
    end
    methods (Access = public)
        function obj = ProjectedFieldComputer(cParams)
            obj.inputData(cParams);
        end
        function compute(obj)
            obj.project();
        end
    end
    methods (Access = private)
        function inputData(obj,cParams)
            obj.beta = cParams.beta;
            obj.eta =cParams.eta;
            obj.filteredField =cParams.filteredField;
        end
    end
    methods (Access = protected,Abstract)
        project(obj) % abstract method
    end
end