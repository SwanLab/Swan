classdef HeavisideProjector < handle

    properties (Access = public)
        projectedField
        derivatedProjectedField
    end

    properties (Access = private)
        beta
        eta
        filteredField
    end

    methods (Access = public)

        function obj = HeavisideProjector(cParams)
            obj.inputData(cParams);
        end

        function updateFilteredField(obj,x_filter)
            obj.filteredField = x_filter;
        end

        function project(obj)
            obj.projectedField = (tanh(obj.beta*obj.eta) + tanh(obj.beta*(obj.filteredField-obj.eta)))/(tanh(obj.beta*obj.eta) + tanh(obj.beta*(1-obj.eta)));
        end

        function derive(obj)
            obj.derivatedProjectedField = obj.beta*(1 - (tanh(obj.beta*(obj.filteredField-obj.eta))).^2)/(tanh(obj.beta*obj.eta) + tanh(obj.beta*(1-obj.eta)));    
        end

    end
    methods (Access = private)

        function inputData(obj,cParams)
            obj.beta = cParams.beta;
            obj.eta  = cParams.eta;
        end
        
    end    
end