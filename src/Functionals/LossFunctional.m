classdef LossFunctional < handle

    properties (Access = private)       
        iBatch                
        order
        nBatches
    end

    properties (Access = private)
        costType
        designVariable
        network
        data
    end
    
    methods (Access = public)
        
        function obj = LossFunctional(cParams)
            obj.init(cParams);
            obj.iBatch = 1;            
            obj.computeNumberOfBatchesAndOrder();
        end

        function [j, dj] = computeFunctionAndGradient(obj, x)
            obj.designVariable.thetavec = x;
            Xb = obj.data.Xtrain;
            Yb = obj.data.Ytrain;
            yOut = obj.network.computeYOut(Xb);
            j  = obj.computeCost(yOut,Yb);
            dj = obj.computeGradient(yOut,Yb);
        end

        function [j,dj,isBD] = computeStochasticCostAndGradient(obj,x,moveBatch)
            obj.designVariable.thetavec = x;            
            Xt = obj.data.Xtrain;
            Yt = obj.data.Ytrain;            
            [Xb,Yb] = obj.updateSampledDataSet(Xt,Yt,obj.iBatch);
            yOut = obj.network.computeYOut(Xb);
            j  = obj.computeCost(yOut,Yb);
            dj = obj.computeGradient(yOut,Yb);
            obj.iBatch = obj.updateBatchCounter(obj.iBatch,moveBatch);
            isBD = obj.isBatchDepleted(obj.iBatch,moveBatch);
        end

        function testError = getTestError(obj)
            Xtest = obj.data.Xtest;
            Ytest = obj.data.Ytest;
            [~,Ypred] = max(obj.network.computeLastH(Xtest), [], 2);
            [~,Ytarget]   = max(Ytest,[],2);
            testError = mean(Ypred ~= Ytarget);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.costType       = cParams.costType;
            obj.designVariable = cParams.designVariable;
            obj.network        = cParams.network;
            obj.data           = cParams.data;
        end

        function j = computeCost(obj,yOut,Yb)
            [j, ~] = obj.lossFunction(Yb,yOut);
        end

        function dj = computeGradient(obj,yOut,Yb)
            [~, dLF] = obj.lossFunction(Yb,yOut);
            dj = obj.network.backprop(Yb,dLF);
        end

        function [J,gc] = lossFunction(obj,y,yOut)
            type = obj.costType;
            yp = yOut - 1e-11;
            switch type
                case '-loglikelihood'
                    c = sum((1-y).*(-log(1-yp)) + y.*(-log(yp)),2);
                    J = mean(c);
                    gc = (yp-y)./(yp.*(1-yp));
                case 'L2'
                    c = ((yp-y).^2);
                    J = sqrt(sum(c, 'all'));
                    gc = (yp-y);
                otherwise
                    msg = [type,' is not a valid cost function'];
                    error(msg)
            end
        end

        function [ord,nBatches] = computeNumberOfBatchesAndOrder(obj)
            nD = size(obj.data.Xtrain,1);            
            [batchSize] = obj.computeBatchSize();
            nBatches = fix(nD/batchSize);
            if nBatches == 1 || nBatches == 0
                ord = 1:nD;
                nBatches = 1;
            else
                ord = randperm(nD,nD);
            end
            obj.order = ord;
            obj.nBatches = nBatches;
        end

        function itIs = isBatchDepleted(obj,iBatch,moveBatch)
            itIs = iBatch == obj.nBatches && moveBatch;
        end

        function iB = updateBatchCounter(obj,iB,moveBatch)
            if iB < obj.nBatches && moveBatch
                iB = iB + 1;
            elseif iB == obj.nBatches && moveBatch
                iB = 1;             
            end
        end

        function batchSize = computeBatchSize(obj)
            if size(obj.data.Xtrain,1) > 200
                batchSize = 200;
            else
                batchSize = size(obj.data.Xtrain,1);
            end
        end        

        function [x,y] = updateSampledDataSet(obj,Xl,Yl,iBatch)
            iB = iBatch;            
            batchSize = obj.computeBatchSize();
            nB = batchSize;
            if iB == fix(size(Xl,1)/nB)
                plus = mod(size(Xl,1),nB);
                x = zeros([nB+plus,size(Xl,2)]);
                y = zeros([nB+plus,size(Yl,2)]);
            else
                plus = 0;
                x = zeros([nB,size(Xl,2)]);
                y = zeros([nB,size(Yl,2)]);
            end
            cont = 1;
            for jB = (iB-1)*nB+1:iB*nB+plus
                x(cont,:) = Xl(obj.order(jB),:);
                y(cont,:) = Yl(obj.order(jB),:);
                cont = cont+1;
            end
        end

    end
end