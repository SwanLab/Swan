classdef OptimizationProblem < handle
    
    properties (Access = private)
        data
        networkParams
        optimizerParams
        costParams
        network
        costFunc
        designVariable
        optimizer
        plotter
        loss
        regularization
    end
  
   methods (Access = public)

       function obj = OptimizationProblem(cParams)
           obj.init(cParams);
           obj.createNetwork();
           obj.createDesignVariable();
           obj.createLossFunctional();
           obj.createRegularizationFunctional();
           obj.createCost();
           obj.createPlotter();
           obj.createOptimizer();
       end

       function solve(obj)
           obj.optimizer.compute();
       end

       function plotCostFnc(obj)
           obj.optimizer.plotCostFunc();
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

       function plotSurface(obj)
           obj.plotter.drawSurfaceResults();
       end

       function plotImage(obj,row)
           obj.plotter.image(row);
       end
       
       function E = computeError(obj,X,Y)
           E = obj.network.forwardprop(X,Y);
       end
       
       function yOut = computeOutputValues(obj,X)
           yOut = obj.network.computeYOut(X);
       end

       function dY = computeGradient(obj,X)
           dY = obj.network.networkGradient(X);
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
           dv = obj.network.getLearnableVariables();
           obj.designVariable = dv;
       end

       function createLossFunctional(obj)
            s.network        = obj.network;
            s.designVariable = obj.designVariable;
            s.data           = obj.data;
            l = Sh_Func_Loss(s);
            obj.loss = l;
        end

        function createRegularizationFunctional(obj)
            s.designVariable = obj.designVariable;
            r = Sh_Func_L2norm(s);
            obj.regularization = r;
        end

       function createCost(obj)
           s.shapeFunctions = {obj.loss, obj.regularization};
           s.weights = [1, obj.costParams.lambda];
           s.Msmooth = [];
           obj.costFunc = Cost(s);
       end

       function createOptimizer(obj)
           s             = obj.optimizerParams;
           s.costFunc    = obj.costFunc;
           s.designVariable = obj.network.getLearnableVariables();
           s.type        = 'SGD';
           s.data        = obj.data;
           s.Xtrain = obj.data.Xtrain;
           s.Ytrain = obj.data.Ytrain;
           s.Xtest = obj.data.Xtest;
           s.Ytest = obj.data.Ytest;
           s.plotter = obj.plotter;
           op = Trainer.create(s);
           obj.optimizer = op;
       end

       function createPlotter(obj)
           s.network   = obj.network;
           s.data      = obj.data;
           s.costfunc  = obj.costFunc;
           obj.plotter = PlotterNN(s);
       end
         
   end

end