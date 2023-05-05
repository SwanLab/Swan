classdef Network < handle
 
    properties (Access = private)
    end

    properties (GetAccess = public, SetAccess = private)
       neuronsPerLayer
       nLayers
       lambda
    end

    properties (Access = private)
       delta
       a_fcn
       HUtype
       OUtype
       Costtype
    end

    properties (Access = private)
       hiddenLayers
       data   
    end

   methods (Access = public)

       function obj = Network(cParams)
           obj.init(cParams);
       end

        function th = computeInitialTheta(obj)
           nPL    = obj.neuronsPerLayer;
           th     = [];
           for i = 2:obj.nLayers
                if i ~= obj.nLayers
                    getB = zeros([1,nPL(i)]) + 0.1;
                else
                    getB = zeros([1,nPL(i)]) + 1/nPL(i);
                end
                u = (6/(nPL(i-1)+nPL(i)))^0.5;
                getW = (unifrnd(-u,u,[1,nPL(i-1)*nPL(i)]));
                th = [th,getW,getB];
           end      
        end

       function c = forwardprop(obj,theta,Xb,Yb)  
           [W,b] = obj.reshapeInLayerForm(theta);                      
           nLy = obj.nLayers;
           a = cell(nLy,1);
           a{1} = Xb;
           for i = 2:nLy
               g_prev = a{i-1};
               Wi = W{i-1};
               bi = b{i-1};
               h = obj.hypothesisfunction(g_prev,Wi,bi);
               [g,~] = obj.actFCN(h,i);
               a{i} = g;
           end
           obj.a_fcn = a;
           [c,~] = obj.costFunction(Yb,a);
       end

       function dc = backprop(obj,theta,Yb)
           [W,b] = obj.reshapeInLayerForm(theta);
           a = obj.a_fcn;
           nPl = obj.neuronsPerLayer;
           nLy = obj.nLayers;
           m = length(Yb);
           deltag = cell(nLy,1);

           dcW = cell(nLy-1,1);           
           dcB = cell(nLy-1,1);
           
          for k = nLy:-1:2    
               [~,a_der] = obj.actFCN(a{k},k); 
               if k == nLy
                   [~,t1] = obj.costFunction(Yb,a);  
                   deltag{k} = t1.*a_der;
               else                    
                   deltag{k} = (W{k}*deltag{k+1}')'.*a_der;
               end
               dcW{k-1} = (1/m)*(a{k-1}'*deltag{k});
               dcB{k-1} = (1/m)*(sum(deltag{k},1));
          end

           dc = [];
           for i = 2:nLy
               aux1 = [reshape(dcW{i-1},[1,nPl(i-1)*nPl(i)]),dcB{i-1}];
               dc  = [dc,aux1];                            
           end
           
       end        

        function g = computeLastH(obj,theta,X)
           nLy = obj.nLayers;
           [W,b] = obj.reshapeInLayerForm(theta);
           h = obj.hypothesisfunction(X,W{1},b{1});
            [g,~] = obj.actFCN(h,2);
            for i = 2:nLy-1
                h = obj.hypothesisfunction(g,W{i},b{i});
                [g,~] = obj.actFCN(h,i+1);
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
                otherwise
                    msg = [type,' is not a valid activation function'];
                    error(msg)
            end
       end
       
   end

   methods (Access = private)

       function init(obj,cParams)
           obj.data         = cParams.data;
           obj.hiddenLayers = cParams.hiddenLayers;
           obj.neuronsPerLayer = [obj.data.nFeatures,obj.hiddenLayers,obj.data.nLabels];
           obj.nLayers = length(obj.neuronsPerLayer);
           if length(cParams) <= 2
               obj.Costtype = '-loglikelihood';
               obj.HUtype = 'ReLU';
               obj.OUtype = 'softmax';
               obj.lambda = 0;
           else
               obj.Costtype = cParams{3};
               obj.HUtype = cParams{4};
               obj.OUtype = cParams{5};
               obj.lambda = cParams{6};
           end
       end


    

       function [J,gc] = costFunction(obj,y,a)
            type = obj.Costtype;
            yp = a{end}-10^-10;
            switch type
                case '-loglikelihood'
                    c = sum((1-y).*(-log(1-yp)) + y.*(-log(yp)),2);
                    J = mean(c);                        
                    gc = (yp-y)./(yp.*(1-yp));
                case 'L2'
                    c = ((yp-y).^2);
                    J = sum(mean(c,1));
                    gc = (yp-y);
                otherwise
                    msg = [type,' is not a valid cost function'];
                    error(msg)
            end
       end
   end 

   methods (Access = private, Static)

       function h = hypothesisfunction(X,W,b)
           h = X*W + b;
       end

       function value = getW(thv,prevL,nextL)
           aux   = thv(1:prevL*nextL);
           value = reshape(aux,[prevL,nextL]);
       end

       function value = getB(thv,prevL,nextL)
           value = thv(prevL*nextL+1:end);
       end

   end

   methods (Access = private)

       function [W,b] = reshapeInLayerForm(obj,thetavec)
            nPL = obj.neuronsPerLayer;
            last = 1;
            b = cell(obj.nLayers-1,1);
            W = cell(obj.nLayers-1,1);
            for i = 2:obj.nLayers
                aux = nPL(i)*nPL(i-1) + nPL(i);
                next = last + aux;
                thetaI = thetavec(last:next-1); 

                prevL = nPL(i-1);
                nextL = nPL(i);

                b{i-1} = obj.getB(thetaI,prevL,nextL);
                W{i-1} = obj.getW(thetaI,prevL,nextL);

                last = next;
            end
       end 

   end

end