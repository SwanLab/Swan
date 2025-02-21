classdef Sh_Func_Loss < handle

    properties(Access = public)

    end
    
    properties (Access = private)
        designVariable
        network
        nBatches
        currentBatch
        moveBatch
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

        function [j, dj] = computeFunctionAndGradient(obj, x)
            obj.designVariable.thetavec = x;
            j  = obj.computeCost();
            dj = obj.computeGradient();
        end

        function batches_depleted = handleStochasticBatch(obj, moveBatch)
            [nD, nB, ~] = obj.fetchBatchSize();
            if nB == 1 || nB == 0
                order = 1:nD;
                nB = 1;
            else
                order = randperm(nD,nD);
            end

            obj.nBatches = nB;
            obj.updateMinibatch(order, obj.currentBatch);
            
            batches_depleted = false;
            if obj.currentBatch < obj.nBatches && moveBatch
                obj.currentBatch = obj.currentBatch + 1;
            elseif obj.currentBatch == obj.nBatches && moveBatch
                obj.currentBatch = 1;
                batches_depleted = true;
            end
        end

        function [nD, nB, batchSize] = fetchBatchSize(obj)
            nD = size(obj.data.Xtrain,1);
            if size(obj.data.Xtrain,1) > 200
                batchSize = 200;
            else
                batchSize = size(obj.data.Xtrain,1);
            end
            nB = fix(nD/batchSize);
        end

        % Funcions getters (eliminar?)
        function [Xtest, Ytest] = getTestData(obj)
            Xtest = obj.data.Xtest;
            Ytest = obj.data.Ytest;
        end

        function h = getOutput(obj,X)
            h = obj.network.computeLastH(X);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.designVariable = cParams.designVariable;
            obj.network        = cParams.network;
            obj.data           = cParams.data;
            obj.currentBatch   = 1;
            obj.moveBatch      = true;
        end

        function j = computeCost(obj)
            j = obj.network.forwardprop(obj.Xb,obj.Yb);
        end

        function dj = computeGradient(obj)
            dj = obj.network.backprop(obj.Yb);
        end
        
        function [x,y,I] = updateMinibatch(obj,order,i)
            Xl = obj.data.Xtrain;
            Yl = obj.data.Ytrain;
            [~, ~, I] = obj.fetchBatchSize();

            cont = 1;
            if i == fix(size(Xl,1)/I)
                plus = mod(size(Xl,1),I);
                x = zeros([I+plus,size(Xl,2)]);
                y = zeros([I+plus,size(Yl,2)]);
            else
                plus = 0;
                x = zeros([I,size(Xl,2)]);
                y = zeros([I,size(Yl,2)]);
            end
            for j = (i-1)*I+1:i*I+plus
                x(cont,:) = Xl(order(j),:);
                y(cont,:) = Yl(order(j),:);
                cont = cont+1;
            end
            obj.Xb = x;
            obj.Yb = y;
        end

    end
end