classdef JuliaNetwork < handle

    properties (GetAccess = public, SetAccess = private)
        neuronsPerLayer
    end

    properties (Access = private)
        struct  % Struct returned from Julia constructor (includes all network internals)
    end

    methods (Access = public)

        function obj = JuliaNetwork(params)
            % Initialize the network by calling Julia constructor
            params.thetavec = [];  % Explicitly pass empty thetavec to avoid double initialization
            obj.struct = callJuliaClass('Network', 'Net', params);

            % Save all input fields in case we need to resend them later
            obj.struct.hiddenLayers    = params.hiddenLayers;
            obj.struct.HUtype          = params.HUtype;
            obj.struct.OUtype          = params.OUtype;
            obj.struct.data            = params.data;
            
        end

        function yOut = computeYOut(obj, Xb)
            params = obj.struct;
            params.Xb = double(Xb');
            result = callJuliaClass('Network', 'computeYOut', params);
            yOut = result.yOut';
        end

        function dc = backprop(obj, Yb, dLF)
            params = obj.struct
            params.Yb = double(Yb);
            params.dLF = double(dLF);
            params.Xb = double(obj.struct.data.Xtrain);
            result = callJuliaClass('Network', 'backprop', params);
            dc = result.dc';
        end

        function dy = networkGradient(obj, X)
            params = obj.struct;
            params.X = double(X);
            result = callJuliaClass('Network', 'networkGradient', params);
            dy = result.dy;
        end

        function H = computeLastH(obj, X)
            params = obj.struct;
            params.X = double(X);
            result = callJuliaClass('Network', 'computeLastH', params);
            H = result.H;
        end

        function vars = getLearnableVariables(obj)
            params = obj.struct;
            vars = callJuliaClass('Network', 'getLearnableVariables', params);
            vars.thetavec = vars.thetavec';  
            obj.struct.thetavec = vars.thetavec; 
            % Struct with fields: thetavec, neuronsPerLayer, nLayers
        end
        
        function setThetavec(obj, x)
        obj.struct.thetavec = x;
        end
        
    end

    

end
