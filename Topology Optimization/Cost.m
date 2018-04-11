classdef Cost < CC
    properties
        weights
    end
    methods
        function obj = Cost(settings,weights,postprocess_TopOpt)
            obj@CC(settings,settings.cost,postprocess_TopOpt);
            if isempty(weights)
                obj.weights = ones(1,length(settings.cost));
            else
                obj.weights = weights;
            end
        end
        
        function updateFields(obj,iSF)
            obj.value = obj.value + obj.weights(iSF)*obj.ShapeFuncs{iSF}.value;
            obj.gradient = obj.gradient + obj.weights(iSF)*obj.ShapeFuncs{iSF}.gradient;
        end
    end
end

