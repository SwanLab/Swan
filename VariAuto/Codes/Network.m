classdef Network < handle
 
    properties (Access = public)
        theta
        cost
        loss
        regularization
        delta
        a_fcn
    end

    properties (GetAccess = public, SetAccess = private)
       neuronsPerLayer
       nLayers
       lambda
       Costtype
       HUtype
       OUtype
    end

    properties (Access = private)
       data  
    end

    properties (Dependent)
        layer
        W
        b
    end

   methods (Access = public)

       function obj = Network(varargin)
           obj.init(varargin);
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

        function [J,gradient] = propagate(obj,theta,Xb,Yb)           
            obj.forwardprop(Xb,Yb);
            J = obj.cost;
            gradient = obj.backprop(theta,Yb);
        end 
       
   end

   methods (Access = private)

       function init(obj,s)
           obj.data = s{1};
           obj.neuronsPerLayer = s{2};
           obj.nLayers = length(s{2});
           if length(s) <= 2
               obj.Costtype = '-loglikelihood';
               obj.HUtype = 'ReLU';
               obj.OUtype = 'softmax';
               obj.lambda = 0;
           else
               obj.Costtype = s{3};
               obj.HUtype = s{4};
               obj.OUtype = s{5};
               obj.lambda = s{6};
           end
       end

       function s = getPrev(obj,thv,prevL,nextL)
            s.theta = thv;
            s.prevL = prevL;
            s.nextL = nextL;
            s.b     = obj.getB(thv,prevL,nextL);
            s.W     = obj.getW(thv,prevL,nextL);
       end

       function forwardprop(obj,Xb,Yb) 
           obj.computeLoss(obj.layer,Xb,Yb);
           obj.computeRegularization(obj.layer);
           c = obj.loss;
           r = obj.regularization;
           l = obj.lambda;
           obj.cost = c + l*r;
       end

       function computeLoss(obj,Xb,Yb)          
           nLy = obj.nLayers;
           a = cell(nLy,1);
           a{1} = Xb;
           for i = 2:nLy
               g_prev = a{i-1};
               obj.W = obj.layer{i-1}.W;
               obj.b = obj.layer{i-1}.b;
               h = obj.hypothesisfunction(g_prev,obj.W,obj.b);
               [g,~] = obj.actFCN(h,i);
               a{i} = g;
           end
           obj.a_fcn = a;
           [J,~] = obj.costFunction(Yb,a);
           obj.loss = J;
       end

       function computeRegularization(obj,layer)
           nLy = obj.nLayers;
           s = 0;
           nth = 0;
           for i = 2:nLy
                s = s + layer{i-1}.theta*layer{i-1}.theta';
                nth = nth + length(layer{i-1}.theta);
           end
           r = 1/(2)*s;
           obj.regularization = r;
       end

       function grad = backprop(obj,theta,Yb)
           obj.layer = obj.getLayer(theta);
           a = obj.a_fcn;
           nPl = obj.neuronsPerLayer;
           nLy = obj.nLayers;
           m = length(Yb);
           deltag = cell(nLy,1);
           gradW = cell(nLy-1,1);
           gradb = cell(nLy-1,1);
           for k = nLy:-1:2    
               [~,a_der] = obj.actFCN(a{k},k); 
               if k == nLy
                   [~,t1] = obj.costFunction(Yb,a);  
                   deltag{k} = t1.*a_der;
               else                    
                   deltag{k} = (obj.layer{k}.W*deltag{k+1}')'.*a_der;
               end
               gradW{k-1} = (1/m)*(a{k-1}'*deltag{k}) + obj.lambda*obj.layer{k-1}.W;
               gradb{k-1} = (1/m)*(sum(deltag{k},1)) + obj.lambda*obj.layer{k-1}.b;
           end
           obj.delta = deltag;
           grad = [];
           for i = 2:nLy
               aux = [reshape(gradW{i-1},[1,nPl(i-1)*nPl(i)]),gradb{i-1}];
               grad = [grad,aux];
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

   methods (Access = private, Static)      
        function h = hypothesisfunction(X,W,b)
          h = X*W + b;
        end     
    end 

   methods 
       function value = getLayer(obj,thetavec)
            nPL = obj.neuronsPerLayer;
            last = 1;
            value = cell(obj.nLayers-1,1);
            for i = 2:obj.nLayers
                aux = nPL(i)*nPL(i-1) + nPL(i);
                next = last + aux;
                theta_i = thetavec(last:next-1);                
                value{i-1} = obj.getPrev(theta_i,nPL(i-1),nPL(i));
                last = next;
            end
       end

       function value = getW(~,thv,prevL,nextL)
            aux = thv(1:prevL*nextL);
            value = reshape(aux,[prevL,nextL]);
       end

        function value = getB(~,thv,prevL,nextL)
            value = thv(prevL*nextL+1:end);
        end
   end

end