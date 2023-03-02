classdef optimizationProblem < handle
 
    properties (GetAccess = public, SetAccess = private)
       data
       thetavec 
       neuronsPerLayer %Network
       nLayers %Network
       lambda %Network
       Costtype
       HUtype %HiddenLayer
       OUtype %OutputLayer
       cost
       regularization
       loss
       gradient
    end

    properties (Dependent)
        layer
    end
    
   properties (Access = private)
       network
       propagator
       plotter
   end

   methods (Access = public)

       function self = optimizationProblem(varargin)
           self.init(varargin);
       end
       
       function computeCost(self)
           self.loss = self.network.loss;
           self.regularization = self.network.regularization;
           self.cost = self.network.cost; 
           self.gradient = self.network.gradient;
       end 

       function plotBoundary(self,type) 
           self.plotter.plotBoundary(type); %No funciona
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
           self.network = Network(self,s); %%%
           self.plotter = Plotter(self);
       end  
   end   
end
