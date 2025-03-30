classdef Network < handle

    properties (GetAccess = public, SetAccess = private)
        neuronsPerLayer
    end

    properties (Access = private)
        delta
        aValues
        HUtype
        OUtype
        Costtype
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

        function c = forwardprop(obj,Xb,Yb)
            obj.computeAvalues(Xb);
            yOut  = obj.aValues{end}; 
            % Out ↓
            [c,~] = obj.lossFunction(Yb,yOut);
            % Out ↑
        end

        function yOut = computeYOut(obj,Xb)
            obj.computeAvalues(Xb);
            yOut = obj.aValues{end}; 
        end

        function computeAvalues(obj,X)
            [W,b] = obj.learnableVariables.reshapeInLayerForm();
            nLy = obj.nLayers;
            z = cell(nLy,1);
            a = cell(nLy,1);
            a{1} = X;
            z{1} = X;
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

        function dc = backprop(obj,Yb)
            [W,b] = obj.learnableVariables.reshapeInLayerForm();
            a = obj.aValues;
            nPl = obj.neuronsPerLayer;
            nLy = obj.nLayers;
            m = length(Yb);
            obj.deltag = cell(nLy,1);

            dcW = cell(nLy-1,1);
            dcB = cell(nLy-1,1);

            for k = nLy:-1:2
                [~,a_der] = obj.actFCN(a{k},k);
                if k == nLy
                    % Out ↓
                    yOut   = a{end};
                    [~,t1] = obj.lossFunction(Yb,yOut);
                    obj.deltag{k} = t1.*a_der;
                    % Out ↑
                else
                    obj.deltag{k} = (W{k}*obj.deltag{k+1}')'.*a_der;
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
            obj.Costtype = cParams.costType;
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
            obj.learnableVariables = t;
        end

        % Out ↓
        function [J,gc] = lossFunction(obj,y,yOut)
            type = obj.Costtype;
            yp = yOut-10^-10;
            switch type
                case '-loglikelihood'
                    c = sum((1-y).*(-log(1-yp)) + y.*(-log(yp)),2);
                    J = mean(c);
                    gc = (yp-y)./(yp.*(1-yp));
                case 'L2'
                    c = ((yp-y).^2);
                    J = sum(mean(c,1));
                    % J = sqrt(sum(c)); % Propuesta 1
                    % J = norm(yp-y,2); % Propuesta 2
                    gc = (yp-y);
                otherwise
                    msg = [type,' is not a valid cost function'];
                    error(msg)
            end
        end
        % Out ↑

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
                    g_der = z.*(1-z);
                case 'ReLU'
                    g = gt(z,0).*z;
                    g_der = gt(z,0);
                case 'tanh'
                    g = (exp(z)-exp(-z))./(exp(z)+exp(-z));
                    g_der = (1-z.^2);
                case 'softmax'
                    g = (exp(z))./(sum(exp(z),2));
                    g_der = z.*(1-z);
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