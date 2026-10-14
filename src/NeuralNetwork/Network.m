classdef Network < handle

    properties (GetAccess = public, SetAccess = private)
        neuronsPerLayer
    end

    properties (Access = private)
        delta
        aValues
        hValues
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

        function [yOut,J] = computeYOutAndJacobian(obj,X)

            obj.computeAvalues(X);

            yOut = obj.aValues{end};

            J = obj.computeJacobian();

        end

        function dc = backprop(obj,Yb,dLF)

            [W,a,h,nLy] = obj.backLoopVars();

            nPl = obj.neuronsPerLayer;

            m = size(Yb,1);

            obj.deltag = cell(nLy,1);

            dcW = cell(nLy-1,1);

            dcB = cell(nLy-1,1);

            for k = nLy:-1:2

                [~,g_der] = obj.actFCN(h{k},k);

                if k == nLy

                    obj.deltag{k} = dLF.*g_der;

                else

                    obj.deltag{k} = ...
                        (W{k}*obj.deltag{k+1}')'.*g_der;

                end

                dcW{k-1} = ...
                    (1/m)*(a{k-1}'*obj.deltag{k});

                dcB{k-1} = ...
                    (1/m)*sum(obj.deltag{k},1);

            end

            dc = [];

            for i = 2:nLy

                aux1 = [...
                    reshape(dcW{i-1},...
                    [1,nPl(i-1)*nPl(i)]),...
                    dcB{i-1}];

                dc = [dc,aux1];

            end

        end

        function J = networkJacobian(obj,X)

            obj.computeAvalues(X);

            J = obj.computeJacobian();

        end

        function g = computeLastH(obj,X)

            nLy = obj.nLayers;

            [W,b] = ...
                obj.learnableVariables.reshapeInLayerForm();

            h = obj.hypothesisfunction(...
                X,W{1},b{1});

            [g,~] = obj.actFCN(h,2);

            for i = 2:nLy-1

                h = obj.hypothesisfunction(...
                    g,W{i},b{i});

                [g,~] = ...
                    obj.actFCN(h,i+1);

            end

        end

        function l = getLearnableVariables(obj)

            l = obj.learnableVariables;

        end

    end

    methods (Access = private)

        function init(obj,cParams)

            obj.hiddenLayers = ...
                cParams.hiddenLayers;

            obj.nFeatures = ...
                cParams.data.nFeatures;

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

            nPL = obj.neuronsPerLayer;

            obj.nLayers = length(nPL);

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

            t = LearnableVariables(s);

            obj.learnableVariables = t;

        end

        function [W,a,h,nLy] = backLoopVars(obj)

            [W,~] = ...
                obj.learnableVariables.reshapeInLayerForm();

            a = obj.aValues;

            h = obj.hValues;

            nLy = obj.nLayers;

        end

        function computeAvalues(obj,X)

            [W,b] = ...
                obj.learnableVariables.reshapeInLayerForm();

            nLy = obj.nLayers;

            a = cell(nLy,1);

            hValues = cell(nLy,1);

            a{1} = X;

            hValues{1} = X;

            for i = 2:nLy

                g_prev = a{i-1};

                Wi = W{i-1};

                bi = b{i-1};

                h = obj.hypothesisfunction(...
                    g_prev,Wi,bi);

                [g,~] = obj.actFCN(h,i);

                hValues{i} = h;

                a{i} = g;

            end

            obj.aValues = a;

            obj.hValues = hValues;

        end

        function J = computeJacobian(obj)

            [W,a,h,nLy] = obj.backLoopVars();

            nSamples = size(a{1},1);

            nFeatures = size(a{1},2);

            J = repmat(...
                eye(nFeatures),...
                1,1,nSamples);

            for k = 2:nLy

                [~,g_der] = ...
                    obj.actFCN(h{k},k);

                J = ...
                    pagemtimes(...
                    J,W{k-1});

                gDer = ...
                    permute(...
                    g_der,[3,2,1]);

                J = J.*gDer;

            end

        end

        function [g,g_der] = actFCN(obj,z,k)

            nLy = obj.nLayers;

            if k == nLy

                type = obj.OUtype;

            else

                type = obj.HUtype;

            end

            switch type

                case 'sigmoid'

                    g = ...
                        1./(1+exp(-z));

                    g_der = ...
                        g.*(1-g);

                case 'ReLU'

                    g = ...
                        gt(z,0).*z;

                    g_der = ...
                        gt(z,0);

                case 'SoftPlus'

                    g = ...
                        max(z,0) + ...
                        log1p(exp(-abs(z)));

                    g_der = ...
                        obj.computeSoftPlusDerivative(z);

                case 'tanh'

                    g = tanh(z);

                    g_der = ...
                        1-g.^2;

                case 'softmax'

                    ez = ...
                        exp(z-max(z,[],2));

                    g = ...
                        ez./sum(ez,2);

                    g_der = ...
                        g.*(1-g);

                case 'linear'

                    g = z;

                    g_der = ...
                        ones(size(z));

                otherwise

                    msg = [...
                        type,...
                        ' is not a valid activation function'];

                    error(msg)

            end

        end

        function g_der = ...
                computeSoftPlusDerivative(~,z)

            g_der = zeros(size(z));

            isPositive = z >= 0;

            g_der(isPositive) = ...
                1./(1+exp(-z(isPositive)));

            ez = exp(z(~isPositive));

            g_der(~isPositive) = ...
                ez./(1+ez);

        end

    end

    methods (Access = private, Static)

        function h = hypothesisfunction(X,W,b)

            h = X*W + b;

        end

    end

end