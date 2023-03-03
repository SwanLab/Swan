classdef FieldProjector < handle
    properties (Access = public)
        projectedField
    end
    properties (Access = private)
        beta 
        eta
        filteredField
        
    end
    methods (Access = public)
        function obj = FieldProjector(cParams)
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
        function project(obj)
            obj.projectedField   = (tanh(obj.beta*obj.eta) + tanh(obj.beta*(obj.filteredField-obj.eta)))/(tanh(obj.beta*obj.eta) + tanh(obj.beta*(1-obj.eta)));
            %volfracD = volumenFraction*sum(xPhysD(:))/sum(xPhysI(:));
        end
    end
end