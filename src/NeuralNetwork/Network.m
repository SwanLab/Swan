classdef Network < handle

    properties (GetAccess = public, SetAccess = private)
        neuronsPerLayer
    end

    properties (Access = private)
        delta
        aValues
        zValues

        HUtype
        OUtype

        nLayers
        hiddenLayers

        learnableVariables

        nFeatures
        nLabels
        nPolyFeatures

        deltag
    end


    methods (Access = public)

        function obj = Network(cParams)

            obj.init(cParams);
            obj.createLearnableVariables();

        end


        function yOut = computeYOut(obj,Xb)

            obj.computeAvalues(Xb);

            yOut = obj.aValues{end};

        end


        function dc = backprop(obj,Yb,dLF)

            %
            % Backpropagation:
            %
            % delta_k =
            %   (delta_{k+1} W_k^T) .* sigma'(z_k)
            %
            % z_k is stored during the forward pass.
            %

            [W,a,z,nLy] = obj.backLoopVars();

            nPl = obj.neuronsPerLayer;

            %
            % Number of samples in the current batch
            %
            m = size(Yb,1);

            obj.deltag = cell(nLy,1);

            dcW = cell(nLy-1,1);
            dcB = cell(nLy-1,1);

            for k = nLy:-1:2

                %
                % Derivative evaluated at the PRE-ACTIVATION z_k
                %
                [~,g_der] = obj.actFCN(z{k},k);

                if k == nLy

                    %
                    % Output layer
                    %
                    obj.deltag{k} = ...
                        dLF .* g_der;

                else

                    %
                    % Hidden layers
                    %
                    obj.deltag{k} = ...
                        (obj.deltag{k+1} * W{k}') ...
                        .* g_der;

                end

                %
                % Gradients with respect to weights and biases.
                %
                % The factor 1/m corresponds to the mean
                % over the samples of the batch.
                %
                dcW{k-1} = ...
                    (a{k-1}' * obj.deltag{k}) / m;

                dcB{k-1} = ...
                    sum(obj.deltag{k},1) / m;

            end


            %
            % Return gradient using the same flattened ordering
            % employed by LearnableVariables.
            %
            dc = [];

            for i = 2:nLy

                aux = [ ...
                    reshape( ...
                        dcW{i-1}, ...
                        [1,nPl(i-1)*nPl(i)]), ...
                    dcB{i-1}];

                dc = [dc,aux]; %#ok<AGROW>

            end

        end


        function dY = networkDirectionalDerivative(obj,X,dX)

            %
            % Forward-mode directional derivative.
            %
            % X:
            %   nPts x nFeatures
            %
            % dX:
            %   nPts x nFeatures x nDir
            %
            % dY:
            %   nPts x nLabels x nDir
            %
            % nDir is completely generic. It may represent
            % derivatives with respect to:
            %
            %   p1, ..., pnVar
            %
            % after the appropriate chain rule has been
            % assembled outside the network.
            %

            nPts = size(X,1);
            nDir = size(dX,3);

            if size(dX,1) ~= nPts

                error( ...
                    'Network:DirectionalDerivative', ...
                    'dX must have the same number of points as X.');

            end

            if size(dX,2) ~= size(X,2)

                error( ...
                    'Network:DirectionalDerivative', ...
                    ['The second dimension of dX must ', ...
                     'coincide with the number of features of X.']);

            end


            [W,b] = ...
                obj.learnableVariables.reshapeInLayerForm();

            nLy = obj.nLayers;

            a  = X;
            da = dX;


            for k = 2:nLy

                Wi = W{k-1};
                bi = b{k-1};

                nIn  = size(Wi,1);
                nOut = size(Wi,2);

                %
                % z_k = a_{k-1} W_k + b_k
                %
                z = ...
                    obj.hypothesisfunction(a,Wi,bi);

                %
                % a_k = sigma(z_k)
                %
                [aNext,g_der] = ...
                    obj.actFCN(z,k);


                %
                % dz_k = da_{k-1} W_k
                %
                % Reshape only combines point/direction so that
                % all points and all directions are processed
                % simultaneously.
                %
                daMat = permute(da,[1 3 2]);

                daMat = ...
                    reshape( ...
                        daMat, ...
                        nPts*nDir, ...
                        nIn);

                dzMat = daMat * Wi;

                dz = ...
                    reshape( ...
                        dzMat, ...
                        nPts, ...
                        nDir, ...
                        nOut);

                dz = permute(dz,[1 3 2]);


                %
                % da_k = sigma'(z_k) .* dz_k
                %
                da = dz .* g_der;

                a = aNext;

            end

            dY = da;

        end


        function J = networkJacobian(obj,X)

            %
            % Full Jacobian of the network output with respect
            % to the NETWORK INPUT FEATURES.
            %
            % Output:
            %
            % J(point,output,inputFeature)
            %

            nPts    = size(X,1);
            nLabels = obj.nLabels;

            %
            % Forward pass stores both a_k and z_k
            %
            obj.computeAvalues(X);

            [W,~,z,nLy] = obj.backLoopVars();


            %
            % Start from d y / d y = I
            %
            J = repmat( ...
                eye(nLabels), ...
                [1,1,nPts]);

            J = permute(J,[3,1,2]);


            %
            % Reverse chain rule
            %
            for k = nLy:-1:2

                %
                % Always use sigma'(z_k), including the
                % output layer.
                %
                % For a linear output this naturally gives 1.
                %
                [~,g_der] = ...
                    obj.actFCN(z{k},k);

                J = ...
                    J .* ...
                    reshape( ...
                        g_der, ...
                        nPts, ...
                        size(g_der,2), ...
                        1);

                J = ...
                    pagemtimes( ...
                        J, ...
                        W{k-1}');

            end


            %
            % nPts x nLabels x nInputFeatures
            %
            J = permute(J,[1,3,2]);

        end


        function g = computeLastH(obj,X)

            obj.computeAvalues(X);

            %
            % Activation of the last hidden layer.
            %
            g = obj.aValues{end-1};

        end


        function l = getLearnableVariables(obj)

            l = obj.learnableVariables;

        end

    end


    methods (Access = private)

        function init(obj,cParams)

            obj.hiddenLayers = ...
                cParams.hiddenLayers;

            %
            % Number of original physical/input variables.
            %
            obj.nFeatures = ...
                cParams.data.nFeatures;

            %
            % Actual number of features entering the network.
            %
            % At the moment this may include polynomial
            % features generated outside Network.
            %
            obj.nPolyFeatures = ...
                size(cParams.data.Xtrain,2);

            obj.nLabels = ...
                cParams.data.nLabels;

            obj.createNeuronsPerLayer();
            obj.createNumberOfLayers();

            obj.HUtype = ...
                cParams.HUtype;

            obj.OUtype = ...
                cParams.OUtype;

        end


        function createNumberOfLayers(obj)

            obj.nLayers = ...
                length(obj.neuronsPerLayer);

        end


        function createNeuronsPerLayer(obj)

            nF = obj.nPolyFeatures;
            hL = obj.hiddenLayers;
            nL = obj.nLabels;

            obj.neuronsPerLayer = ...
                [nF,hL,nL];

        end


        function createLearnableVariables(obj)

            s.neuronsPerLayer = ...
                obj.neuronsPerLayer;

            s.nLayers = ...
                obj.nLayers;

            obj.learnableVariables = ...
                LearnableVariables(s);

        end


        function [W,a,z,nLy] = backLoopVars(obj)

            [W,~] = ...
                obj.learnableVariables.reshapeInLayerForm();

            a = obj.aValues;
            z = obj.zValues;

            nLy = obj.nLayers;

        end


        function computeAvalues(obj,X)

            %
            % Forward propagation.
            %
            % We deliberately store BOTH:
            %
            %   z_k = pre-activation
            %   a_k = activation
            %
            % because the mathematical derivative of a layer
            % is naturally sigma'(z_k).
            %

            [W,b] = ...
                obj.learnableVariables.reshapeInLayerForm();

            nLy = obj.nLayers;

            a = cell(nLy,1);
            z = cell(nLy,1);

            %
            % Input layer
            %
            a{1} = X;
            z{1} = [];


            for k = 2:nLy

                z{k} = ...
                    obj.hypothesisfunction( ...
                        a{k-1}, ...
                        W{k-1}, ...
                        b{k-1});

                a{k} = ...
                    obj.actFCN(z{k},k);

            end

            obj.aValues = a;
            obj.zValues = z;

        end


        function [g,g_der] = actFCN(obj,z,k)

            %
            % Activation and its derivative with respect
            % to the PRE-ACTIVATION z.
            %

            if k == obj.nLayers

                type = obj.OUtype;

            else

                type = obj.HUtype;

            end


            switch type

                case 'sigmoid'

                    g = ...
                        1 ./ (1 + exp(-z));

                    if nargout > 1

                        g_der = ...
                            g .* (1-g);

                    end


                case 'ReLU'

                    g = max(z,0);

                    if nargout > 1

                        %
                        % ReLU is not differentiable at z = 0.
                        % We use the standard convention
                        % sigma'(0) = 0.
                        %
                        g_der = ...
                            double(z > 0);

                    end


                case 'tanh'

                    g = tanh(z);

                    if nargout > 1

                        g_der = ...
                            1 - g.^2;

                    end


                case 'linear'

                    g = z;

                    if nargout > 1

                        g_der = ...
                            ones(size(z));

                    end


                case 'softmax'

                    %
                    % Forward softmax is allowed.
                    %
                    % Its derivative, however, is NOT an
                    % element-wise quantity:
                    %
                    % ds_i/dz_j =
                    % s_i (delta_ij - s_j)
                    %
                    % Therefore the current generic element-wise
                    % derivative machinery cannot represent it.
                    %

                    zShift = ...
                        z - max(z,[],2);

                    expZ = exp(zShift);

                    g = ...
                        expZ ./ sum(expZ,2);

                    if nargout > 1

                        error( ...
                            'Network:SoftmaxDerivative', ...
                            ['The softmax derivative is a full ', ...
                             'Jacobian and cannot be represented ', ...
                             'by the current element-wise ', ...
                             'activation derivative interface.']);

                    end


                otherwise

                    error( ...
                        'Network:InvalidActivation', ...
                        '%s is not a valid activation function.', ...
                        type);

            end

        end

    end


    methods (Access = private, Static)

        function h = hypothesisfunction(X,W,b)

            h = X*W + b;

        end

    end

end