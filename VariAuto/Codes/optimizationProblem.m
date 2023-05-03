classdef optimizationProblem < handle
 
    properties (GetAccess = public, SetAccess = private)
       data
       network
       costfnc
    end
    
    properties (Access = private)
       designVariable
       structure
       gradient
       cost
       regularization
       loss
       learningRate
       plotter
       lambda
       delta
       actFCN
    end

   methods (Access = public)

       function obj = optimizationProblem(varargin)
           obj.init(varargin);
           obj.createNetwork();
           obj.createDesignVariable();
           obj.createCost();
           obj.createPlotter();
           optimizer = Trainer.create(obj.costfnc,'SGD',obj.learningRate,obj.data);
           optimizer.train();
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
           s.designVariable = obj.designVariable;
           obj.costfnc = CostFunction(s);
       end

       function createPlotter(obj)
           obj.plotter = Plotter(obj);
       end
         
   end

end