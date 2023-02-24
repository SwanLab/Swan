classdef Propagator < handle
    
    properties (Access = public)
        lambda
    end
    
    properties (GetAccess = public, SetAccess = private)
        regularization
        loss
        cost
    end

    properties (Access = private)
        data
        network
        a_fcn
        delta
    end
    
    methods (Access = public)

        function self = Propagator(data,lambda,net)
            self.data = data;
            self.lambda = lambda;
            self.network = net;
        end

        function [J,gradient] = propagate(self,layer,Xb,Yb)
            self.forwardprop(layer,Xb,Yb);
            J = self.cost;
            gradient = self.backprop(layer,Yb);
        end       

       function g = compute_last_H(self,X)
           nLy = self.network.nLayers;
           layer = self.network.layer;
           h = self.hypothesisfunction(X,layer{1}.W,layer{1}.b);
            [g,~] = self.actFCN(h,2);
            for i = 2:nLy-1
                h = self.hypothesisfunction(g,layer{i}.W,layer{i}.b);
                [g,~] = self.actFCN(h,i+1);
            end
       end      
    end

    methods (Access = private)

       function forwardprop(self,layer,Xb,Yb)
           self.computeLoss(layer,Xb,Yb);
           self.computeRegularization(layer);
           c = self.loss;
           r = self.regularization;
           l = self.lambda;
           self.cost = c + l*r;
       end

       function computeLoss(self,layer,Xb,Yb)          
           nLy = self.network.nLayers;
           a = cell(nLy,1);
           a{1} = Xb;
           for i = 2:nLy
               g_prev = a{i-1};
               h = self.hypothesisfunction(g_prev,layer{i-1}.W,layer{i-1}.b);
               [g,~] = self.actFCN(h,i);
               a{i} = g;
           end
           self.a_fcn = a;
           [J,~] = self.costFunction(Yb,a);
           self.loss = J;
       end

       function computeRegularization(self,layer)
           nLy = self.network.nLayers;
           s = 0;
           nth = 0;
           for i = 2:nLy
                s = s + layer{i-1}.theta*layer{i-1}.theta';
                nth = nth + length(layer{i-1}.theta);
           end
           r = 1/(2)*s;
           self.regularization = r;
       end

       function grad = backprop(self,layer,Yb)
           a = self.a_fcn;
           nPl = self.network.neuronsPerLayer;
           nLy = self.network.nLayers;
           m = length(Yb);
           deltag = cell(nLy,1);
           gradW = cell(nLy-1,1);
           gradb = cell(nLy-1,1);
           for k = nLy:-1:2    
               [~,a_der] = self.actFCN(a{k},k); 
               if k == nLy
                   [~,t1] = self.costFunction(Yb,a);  
                   deltag{k} = t1.*a_der;
               else                    
                   deltag{k} = (layer{k}.W*deltag{k+1}')'.*a_der;
               end
               gradW{k-1} = (1/m)*(a{k-1}'*deltag{k}) + self.lambda*layer{k-1}.W;
               gradb{k-1} = (1/m)*(sum(deltag{k},1)) + self.lambda*layer{k-1}.b;
           end
           self.delta = deltag;
           grad = [];
           for i = 2:nLy
               aux = [reshape(gradW{i-1},[1,nPl(i-1)*nPl(i)]),gradb{i-1}];
               grad = [grad,aux];
           end
       end

       function [J,gc] = costFunction(self,y,a)
            type = self.network.Costtype;
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

       function [g,g_der] = actFCN(self,z,k)
            nLy = self.network.nLayers;
            if k == nLy
                type = self.network.OUtype;
            else
                type = self.network.HUtype;
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
end