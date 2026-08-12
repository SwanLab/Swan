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
            [W, a, nLy] = obj.backLoopVars();
            nPl = obj.neuronsPerLayer;
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

        function dY = networkDirectionalDerivative(obj,X,dX)
            nPts = size(X,1);
            nDir = size(dX,3);           
            if size(dX,1) ~= nPts
                error('dX must have the same number of points as X.');
            end
            if size(dX,2) ~= size(X,2)
                error('The second dimension of dX must coincide with the features of X.');
            end           
            [W,b] = obj.learnableVariables.reshapeInLayerForm();
            nLy = obj.nLayers;
            a  = X;
            da = dX;
            for k = 2:nLy
                Wi = W{k-1};
                bi = b{k-1};
                nIn  = size(Wi,1);
                nOut = size(Wi,2);
                h = obj.hypothesisfunction(a,Wi,bi);
                [aNext,g_der] = obj.actFCN(h,k);
                daMat = permute(da,[1 3 2]);
                daMat = reshape(daMat,nPts*nDir,nIn);
                dhMat = daMat*Wi;
                dh = reshape(dhMat,nPts,nDir,nOut);
                dh = permute(dh,[1 3 2]);
                da = dh .* g_der;
                a = aNext;
            end
            dY = da;
        end
       
        function J = networkJacobian(obj, X)
            nPts = size(X, 1);
            nFeatures = obj.nFeatures;
            nLabels = obj.nLabels;
            obj.computeAvalues(X);
            [W, a, nLy] = obj.backLoopVars();           
            J = repmat(eye(nLabels), [1, 1, nPts]);
            J = permute(J, [3, 1, 2]);  
            for k = nLy:-1:2
                if k == nLy
                    g_der = ones(size(a{k}));  
                else
                    [~, g_der] = obj.actFCN(a{k}, k);     
                end                
                J = J .* reshape(g_der, nPts, size(g_der, 2), 1);                
                J = pagemtimes(J, W{k-1}');
            end            
            J = permute(J, [1, 3, 2]);
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
            obj.learnableVariables = t;
        end

        function [W, a, nLy] = backLoopVars(obj)
            [W,~] = obj.learnableVariables.reshapeInLayerForm();
            a     = obj.aValues;
            nLy   = obj.nLayers;
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
                    g_der = z.*(1-z);
                case 'ReLU'
                    g = gt(z,0).*z;
                    g_der = gt(z,0);
                case 'tanh'
                    g     = tanh(z);
                    g_der = 1 - g.^2;
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