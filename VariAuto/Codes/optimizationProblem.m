classdef optimizationProblem < handle
 
    properties (GetAccess = public, SetAccess = private)
       data
       thetavec
       cost
       regularization
       loss
       gradient
       structure
       network
       designVariable
       costfnc
       learningRate
    end
    
    properties (Access = private)
       plotter
       lambda
       a_fcn
       delta
    end

   methods (Access = public)

       function obj = optimizationProblem(varargin)
           obj.init(varargin);
           obj.createNetwork();
           obj.createDesignVariable();
           obj.createPlotter();
           obj.thetavec = obj.designVariable.thetavec;
%           obj.costfnc = costFnc(obj.network);
           optimizer = Trainer.create(obj,'SGD',obj.learningRate);
           optimizer.train();
       end
       
       function computeCost(obj,theta,Xb,Yb)
          obj.network.thetavec = theta;
      %     obj.thetavec = theta;
           [J,grad] = obj.propagate(obj.network.getLayer(),Xb,Yb); 
           obj.loss = obj.loss;
           l = obj.lambda;
           obj.regularization = l*obj.regularization;
           obj.cost = J; 
           obj.gradient = grad;
       end 
              
       function h = getOutput(obj,X)
            h = obj.compute_last_H(X);
       end

       function plotBoundary(obj,type) 
           obj.plotter.plotBoundary(type);
       end

       function plotConections(obj)
           obj.plotter.plotNetworkStatus();
       end

       function plotConfusionMatrix(obj)
           obj.plotter.drawConfusionMat();
       end
   end

   methods (Access = private)

       function init(obj,s)
           obj.data = s{1};
           obj.structure = s{2};
           obj.learningRate = s{3}; 
           obj.lambda = s{4};
       end  

       function createNetwork(obj)
           s1 = obj.data;
           s2 = obj.structure;
           n  = Network(s1,s2);
           obj.network = n;
       end

       function createDesignVariable(obj)
           s.initValue = obj.network.computeInitialTheta();
           t = DesVar(s);
           obj.designVariable = t;
       end

       function createPlotter(obj)
           obj.plotter = Plotter(obj);
       end

       function [J,gradient] = propagate(obj,layer,Xb,Yb)
            obj.forwardprop(layer,Xb,Yb);
            J = obj.cost;
            gradient = obj.backprop(layer,Yb);
       end  

       function g = compute_last_H(obj,X)
           nLy = obj.network.nLayers;
           layer = obj.network.getLayer();
           h = obj.hypothesisfunction(X,layer{1}.W,layer{1}.b);
            [g,~] = obj.actFCN(h,2);
            for i = 2:nLy-1
                h = obj.hypothesisfunction(g,layer{i}.W,layer{i}.b);
                [g,~] = obj.actFCN(h,i+1);
            end
       end   

       function forwardprop(obj,layer,Xb,Yb)
           obj.computeLoss(layer,Xb,Yb);
           obj.computeRegularization(layer);
           c = obj.loss;
           r = obj.regularization;
           l = obj.lambda;
           obj.cost = c + l*r;
       end

       function computeLoss(obj,layer,Xb,Yb)          
           nLy = obj.network.nLayers;
           a = cell(nLy,1);
           a{1} = Xb;
           for i = 2:nLy
               g_prev = a{i-1};
               h = obj.hypothesisfunction(g_prev,layer{i-1}.W,layer{i-1}.b);
               [g,~] = obj.actFCN(h,i);
               a{i} = g;
           end
           obj.a_fcn = a;
           [J,~] = obj.costFunction(Yb,a);
           obj.loss = J;
       end

       function computeRegularization(obj,layer)
           nLy = obj.network.nLayers;
           s = 0;
           nth = 0;
           for i = 2:nLy
                s = s + layer{i-1}.theta*layer{i-1}.theta';
                nth = nth + length(layer{i-1}.theta);
           end
           r = 1/(2)*s;
           obj.regularization = r;
       end

       function grad = backprop(obj,layer,Yb)
           a = obj.a_fcn;
           nPl = obj.network.neuronsPerLayer;
           nLy = obj.network.nLayers;
           m = length(Yb);
           deltag = cell(nLy,1);
           gradW = cell(nLy-1,1);
           gradb = cell(nLy-1,1);
           l = obj.lambda;
           for k = nLy:-1:2    
               [~,a_der] = obj.actFCN(a{k},k); 
               if k == nLy
                   [~,t1] = obj.costFunction(Yb,a);  
                   deltag{k} = t1.*a_der;
               else                    
                   deltag{k} = (layer{k}.W*deltag{k+1}')'.*a_der;
               end
               gradW{k-1} = (1/m)*(a{k-1}'*deltag{k}) + l*layer{k-1}.W;
               gradb{k-1} = (1/m)*(sum(deltag{k},1)) + l*layer{k-1}.b;
           end
           obj.delta = deltag;
           grad = [];
           for i = 2:nLy
               aux = [reshape(gradW{i-1},[1,nPl(i-1)*nPl(i)]),gradb{i-1}];
               grad = [grad,aux];
           end
       end

       function [J,gc] = costFunction(obj,y,a)
            type = obj.network.Costtype;
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
            nLy = obj.network.nLayers;
            if k == nLy
                type = obj.network.OUtype;
            else
                type = obj.network.HUtype;
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