classdef Network < handle

    properties (GetAccess = public, SetAccess = private)
        neuronsPerLayer
    end

    properties (Access = private)
        delta
        aValues
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
            [W,~] = obj.learnableVariables.reshapeInLayerForm();
            a = obj.aValues;
            nPl = obj.neuronsPerLayer;
            nLy = obj.nLayers;
            m = length(Yb);
            obj.deltag = cell(nLy,1);
            dcW = cell(nLy-1,1);
            dcB = cell(nLy-1,1);
            for k = nLy:-1:2
                [~,g_der] = obj.actFCN(a{k},k);
                if k == nLy
                    obj.deltag{k} = dLF.*g_der;
                else
                    obj.deltag{k} = (W{k}*obj.deltag{k+1}')'.*g_der;
                end
                dcW{k-1} = (1/m)*(a{k-1}'*obj.deltag{k});
                dcB{k-1} = (1/m)*(sum(obj.deltag{k},1));
            end
            dc = [];
            for i = 2:nLy
                aux1 = [reshape(dcW{i-1},[1,nPl(i-1)*nPl(i)]),dcB{i-1}];
                dc   = [dc,aux1];
            end
        end

        function dy = networkGradient(obj,X)
            obj.computeAvalues(X);
            [W,~] = obj.learnableVariables.reshapeInLayerForm();
            a = obj.aValues;
            nLy = obj.nLayers;
            for k = nLy-1:-1:1
                [~,a_der] = obj.actFCN(a{k+1},k+1);
                a_der = diag(a_der);
                parDer = a_der * W{k}';
                if k == nLy-1
                    grad = parDer;
                else
                    grad = grad * parDer;
                end
            end
            dy = grad;
        end

        function g = computeLastH(obj,X)
            nLy = obj.nLayers;
            [W,b] = obj.learnableVariables.reshapeInLayerForm();
            h = obj.hypothesisfunction(X,W{1},b{1});
            [g,~] = obj.actFCN(h,2);
            for i = 2:nLy-1
                h = obj.hypothesisfunction(g,W{i},b{i});
                [g,~] = obj.actFCN(h,i+1);
            end
        end

        function l = getLearnableVariables(obj)
            l = obj.learnableVariables;
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.hiddenLayers = cParams.hiddenLayers;
            obj.nFeatures    = cParams.data.nFeatures;
            obj.nPolyFeatures= size(cParams.data.Xtrain,2);
            obj.nLabels      = cParams.data.nLabels;
            obj.createNeuronsPerLayer();
            obj.createNumberOfLayers();
            obj.HUtype   = cParams.HUtype;
            obj.OUtype   = cParams.OUtype;
        end

        function createNumberOfLayers(obj)
            nPL = obj.neuronsPerLayer;
            obj.nLayers = length(nPL);
        end
        
        function createNeuronsPerLayer(obj)
            nF = obj.nPolyFeatures;
            hL = obj.hiddenLayers;
            nL = obj.nLabels;
            obj.neuronsPerLayer = [nF,hL,nL];
        end

        function createLearnableVariables(obj)
            s.neuronsPerLayer = obj.neuronsPerLayer;
            s.nLayers         = obj.nLayers;
            t = LearnableVariables(s);
            %t = JuliaLearnableVariables(s);
            obj.learnableVariables = t;
        end

        function computeAvalues(obj,X)
            [W,b] = obj.learnableVariables.reshapeInLayerForm();
            nLy = obj.nLayers;
            a = cell(nLy,1);
            a{1} = X;
            for i = 2:nLy
                g_prev = a{i-1};
                Wi = W{i-1};
                bi = b{i-1};
                h = obj.hypothesisfunction(g_prev,Wi,bi);
                [g,~] = obj.actFCN(h,i);
                a{i} = g;
            end
            obj.aValues = a;
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
                    g = 1./(1+exp(-z));
                    g_der = g.*(1-g);
                case 'ReLU'
                    g = gt(z,0).*z;
                    g_der = gt(z,0);
                case 'tanh'
                    g = (exp(z)-exp(-z))./(exp(z)+exp(-z));
                    g_der = (1-g.^2);
                case 'softmax'
                    g = (exp(z))./(sum(exp(z),2));
                    g_der = g.*(1-g);
                case 'linear'
                    g = z;
                    g_der = ones(size(z));
                otherwise
                    msg = [type,' is not a valid activation function'];
                    error(msg)
            end
        end

    end

    methods (Access = private, Static)

        function h = hypothesisfunction(X,W,b)
            h = X*W + b;
        end

    end

end