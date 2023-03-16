classdef FieldProjector < ProjectedFieldComputer
    properties (Access = public)
        projectedField
    end
    methods (Access = protected)
        function project(obj)
            obj.projectedField   = (tanh(obj.beta*obj.eta) + tanh(obj.beta*(obj.filteredField-obj.eta)))/(tanh(obj.beta*obj.eta) + tanh(obj.beta*(1-obj.eta)));
            %volfracD = volumenFraction*sum(xPhysD(:))/sum(xPhysI(:));
        end
    end
end