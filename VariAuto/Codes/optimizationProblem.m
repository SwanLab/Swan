classdef optimizationProblem < handle
 
    properties (GetAccess = public, SetAccess = private)
       data
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
           obj.createCost();
           obj.createPlotter();
           optimizer = Trainer.create(obj.costfnc,'SGD',obj.learningRate);
           optimizer.train();
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
           t = DesignVariable(s);
           obj.designVariable = t;
       end

       function createCost(obj)
           s.data = obj.data;
           s.network = obj.network;
           s.lambda  = obj.lambda;
           obj.costfnc = CostFunction(s);
       end


       function createPlotter(obj)
           obj.plotter = Plotter(obj);
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
   end
    
    methods (Access = private, Static)      
        function h = hypothesisfunction(X,W,b)
          h = X*W + b;
        end     
    end  
end