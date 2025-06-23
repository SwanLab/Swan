classdef JuliaNetwork < handle

    properties (GetAccess = public, SetAccess = private)
        neuronsPerLayer
    end

    properties (Access = private)
        data  % Struct returned from Julia constructor (includes all network internals)
    end

    methods (Access = public)

        function obj = JuliaNetwork(params)
            % Initialize the network by calling Julia constructor
            obj.data = callJuliaClass('Network', 'Net', params);

            % Save all input fields in case we need to resend them later
            obj.data.hiddenLayers    = params.hiddenLayers;
            obj.data.nFeatures       = params.data.nFeatures;
            obj.data.nPolyFeatures   = size(params.data.Xtrain,2);
            obj.data.nLabels         = params.data.nLabels;
            obj.data.HUtype          = params.HUtype;
            obj.data.OUtype          = params.OUtype;
        end

        function yOut = computeYOut(obj, Xb)
            params = obj.data;
            params.Xb = double(Xb);
            result = callJuliaClass('Network', 'computeYOut', params);
            yOut = result.yOut;
        end

        function dc = backprop(obj, Yb, dLF)
            params = obj.data;
            params.Yb = double(Yb);
            params.dLF = double(dLF);
            result = callJuliaClass('Network', 'backprop', params);
            dc = result.dc';
        end

        function dy = networkGradient(obj, X)
            params = obj.data;
            params.X = double(X);
            result = callJuliaClass('Network', 'networkGradient', params);
            dy = result.dy;
        end

        function H = computeLastH(obj, X)
            params = obj.data;
            params.X = double(X);
            result = callJuliaClass('Network', 'computeLastH', params);
            H = result.H;
        end

        function vars = getLearnableVariables(obj)
            params = obj.data;
            result = callJuliaClass('Network', 'getLearnableVariables', params);
            vars = result;  % Struct with fields: thetavec, neuronsPerLayer, nLayers
        end

    end
end
