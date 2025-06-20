classdef JuliaLearnableVariables < handle
    properties (Access = private)
        data  % Struct returned from Julia constructor (holds neuronsPerLayer, nLayers, thetavec)
    end
    
    methods (Access = public)
        function obj = JuliaLearnableVariables(params)
            % params should contain 'neuronsPerLayer' and 'nLayers'
            % Call Julia constructor and save returned data (thetavec)
            obj.data = callJuliaClass('LearnableVariables', 'LearnableVars', params);
            
            % Save params fields too, to pass in future calls
            obj.data.neuronsPerLayer = params.neuronsPerLayer;
            obj.data.nLayers = params.nLayers;
        end
        
        function [W, b] = reshapeInLayerForm(obj)
            % Prepare full parameters for the method call
            params.neuronsPerLayer = obj.data.neuronsPerLayer;
            params.nLayers = obj.data.nLayers;
            params.thetavec = obj.data.thetavec;
            
            % Call Julia method
            output = callJuliaClass('LearnableVariables', 'reshapeInLayerForm', params);
            
            % Return outputs (W and b)
            W = output.W;
            b = output.b;
        end
        
    end
end
