classdef optimizationProblem < handle
 
    properties (GetAccess = public, SetAccess = private)
       data
       thetavec
       cost
       regularization
       loss
       gradient
    end
    
   properties (Access = private)
       network
       propagator
       plotter
   end

   methods (Access = public)

       function self = optimizationProblem(varargin)
           self.init(varargin);
           self.thetavec = self.network.thetavec;
       end
       
       function computeCost(self,theta,Xb,Yb)
           self.thetavec = theta;
           [J,grad] = self.propagator.propagate(self.network.layer,Xb,Yb); 
           self.loss = self.propagator.loss;
           l = self.network.lambda;
           self.regularization = l*self.propagator.regularization;
           self.cost = J; 
           self.gradient = grad;
       end 
              
       function h = getOutput(self,X)
            h = self.propagator.compute_last_H(X);
       end

       function plotBoundary(self,type) 
           self.plotter.plotBoundary(type);
       end

       function plotConections(self)
           self.plotter.plotNetworkStatus();
       end

       function plotConfusionMatrix(self)
           self.plotter.drawConfusionMat();
       end
   end

   methods (Access = private)

       function init(self,s)
           self.data = s{1};
           self.network = s{2};
           self.propagator = Propagator(self.network.data,self.network.lambda,self);
           %self.plotter = Plotter(self);
       end  
   end   
end
