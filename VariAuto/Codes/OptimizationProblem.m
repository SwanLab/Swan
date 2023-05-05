classdef OptimizationProblem < handle
 
   properties (Access = private)
       network
       costFunc
       designVariable
       optimizer
       plotter       
    end

    properties (Access = private)
        data
        networkParams
        optimizerParams
        costParams
    end

   methods (Access = public)

       function obj = OptimizationProblem(cParams)
           obj.init(cParams);
           obj.createNetwork();
           obj.createDesignVariable();
           obj.createCost();
           obj.createPlotter();
           obj.createOptimizer();           
       end

       function solve(obj)
           obj.optimizer.train();
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

       function init(obj,cParams)
           obj.data            = cParams.data;
           obj.networkParams   = cParams.networkParams;
           obj.optimizerParams = cParams.optimizerParams; 
           obj.costParams      = cParams.costParams;
       end  

       function createNetwork(obj)
           s      = obj.networkParams;
           s.data = obj.data;
           n  = Network(s);
           obj.network = n;
       end

       function createDesignVariable(obj)
           s.initValue = obj.network.computeInitialTheta();
           t = DesignVariable(s);
           obj.designVariable = t;
       end

       function createCost(obj)
           s         = obj.costParams;
           s.network = obj.network;
           s.designVariable = obj.designVariable;
           obj.costFunc = CostFunction(s);
       end

       function createOptimizer(obj)
           s             = obj.optimizerParams;
           s.costFunc    = obj.costFunc;
           s.designVariable = obj.designVariable;
           s.type        = 'SGD';
           s.data        = obj.data;
           s.maxFunEvals = 2000;
           op = Trainer.create(s);
           obj.optimizer = op;
       end

       function createPlotter(obj)
           s.network   = obj.network;
           s.data      = obj.data;
           s.costfunc  = obj.costFunc;
           obj.plotter = Plotter(s);
       end
         
   end

end