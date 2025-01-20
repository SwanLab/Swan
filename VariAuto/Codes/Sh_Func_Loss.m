classdef Sh_Func_Loss < handle
    
    properties (Access = private)
        designVariable
        network

        % Pau addition
        data
        Xb
        Yb
    end
    
    methods (Access = public)
        
        function obj = Sh_Func_Loss(cParams)
            obj.init(cParams)            
        end
        
        function [j,dj] = computeCostAndGradient(obj)                        
            j  = obj.computeCost();
            dj = obj.computeGradient();            
        end

        function [nD, nB, batchSize] = getBatchSize(obj)
            nD = obj.data.Batch_nD;
            nB = obj.data.Batch_nB;
            batchSize = obj.data.batchSize;
        end

        function [Xtest, Ytest] = getTestData(obj)
            Xtest = obj.data.Xtest;
            Ytest = obj.data.Ytest;
        end

        function obj = updateBatch(obj, order, i)
            [X, Y] = obj.data.createMinibatch(order, i);
            obj.Xb = X;
            obj.Yb = Y;
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.designVariable = cParams.designVariable;
            obj.network        = cParams.network;
            obj.data           = cParams.data;
        end

       function j = computeCost(obj)
           j = obj.network.forwardprop(obj.Xb,obj.Yb);
       end

       function dj = computeGradient(obj)
           dj = obj.network.backprop(obj.Yb);
       end                  
        
    end
    
end