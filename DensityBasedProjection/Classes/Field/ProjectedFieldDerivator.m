classdef ProjectedFieldDerivator < FieldComputer
    properties (Access = public)
        derivatedProjectedField
    end
    methods (Access = protected)
        function project(obj)
            obj.derivatedProjectedField     = obj.beta*(1 - (tanh(obj.beta*(obj.filteredField-obj.eta))).^2)/(tanh(obj.beta*obj.eta) + tanh(obj.beta*(1-obj.eta)));    
        end
    end
end