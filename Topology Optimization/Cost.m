classdef Cost < CC
    
    properties (Access = private)
        weights
    end
    
    methods (Access = public)
        
        function obj = Cost(settings,designVar)
            obj.init(settings, settings.cost,designVar);            
            obj.setWeights(settings.weights);
        end
        
    end
    
    methods (Access = public)
        
        function updateFields(obj,iSF)
            obj.value = obj.value + obj.weights(iSF)*obj.ShapeFuncs{iSF}.value;
            obj.gradient = obj.gradient + obj.weights(iSF)*obj.ShapeFuncs{iSF}.gradient;
        end
        
    end
    
    methods (Access = private)
        
        function setWeights(obj,weights)
            if isempty(weights)
                obj.weights = ones(1,length(settings.cost));
            else
                obj.weights = weights;
            end
        end
        
    end
    
end
