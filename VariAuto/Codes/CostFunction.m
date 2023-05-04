classdef CostFunction < handle
    
    properties (Access = public)
        cost
        designVariable
    end
    
    properties (Access = private)
        data
        network
        lambda
        gradient
        regularization
        loss
    end

    methods (Access = public)

        function obj = CostFunction(cParams)
            obj.init(cParams);
        end
        
        function [J,grad] = computeCost(obj,theta,Xb,Yb)
           obj.designVariable.thetavec = theta;
           [J,grad] = obj.network.propagate(theta,Xb,Yb); 
           obj.loss = obj.network.loss;
           l = obj.network.lambda;
           obj.regularization = l*obj.network.regularization;
           obj.cost = J; 
           obj.gradient = grad;
        end 

        function h = getOutput(obj,X)
            h = obj.compute_last_H(X);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
           obj.data = cParams.data;
           obj.network = cParams.network;
           obj.lambda  = cParams.lambda;
           obj.designVariable = cParams.designVariable;
        end 
        function g = compute_last_H(obj,X)
           nLy = obj.network.nLayers;
           layer = obj.network.getLayer(obj.designVariable.thetavec);
           h = obj.network.hypothesisfunction(X,layer{1}.W,layer{1}.b);
            [g,~] = obj.network.actFCN(h,2);
            for i = 2:nLy-1
                h = obj.network.hypothesisfunction(g,layer{i}.W,layer{i}.b);
                [g,~] = obj.network.actFCN(h,i+1);
            end
        end
       
    end   

end

