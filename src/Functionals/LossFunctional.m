classdef LossFunctional < handle

    properties (Access = private)
        iBatch
        order
        nBatches
        currentBatchIndices
        currentBatchDepleted
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

            obj.currentBatchIndices = [];

            obj.currentBatchDepleted = false;

            obj.computeNumberOfBatchesAndOrder();

        end

        function [j,dj] = computeFunctionAndGradient(obj,x)

            obj.designVariable.thetavec = x;

            Xb = obj.data.Xtrain;
            Yb = obj.data.Ytrain;

            yOut = ...
                obj.network.computeYOut(Xb);

            j = ...
                obj.computeCost(yOut,Yb);

            dj = ...
                obj.computeGradient(yOut,Yb);

        end
        function title = getTitleToPlot(obj)

            title = ...
                ['Loss ',...
                obj.costType];

        end

        function [j,dj] = computeCostAndGradient(obj,x)

            [j,dj] = ...
                obj.computeFunctionAndGradient(x);

        end

        function [j,dj,isBD] = ...
                computeStochasticCostAndGradient(...
                obj,x,moveBatch)

            obj.designVariable.thetavec = x;

            Xt = ...
                obj.data.Xtrain;

            Yt = ...
                obj.data.Ytrain;

            obj.updateCurrentBatch(moveBatch);

            indices = ...
                obj.currentBatchIndices;

            Xb = ...
                Xt(indices,:);

            Yb = ...
                Yt(indices,:);

            yOut = ...
                obj.network.computeYOut(Xb);

            j = ...
                obj.computeCost(yOut,Yb);

            dj = ...
                obj.computeGradient(yOut,Yb);

            isBD = ...
                obj.currentBatchDepleted;

        end

        function testError = getTestError(obj)

            Xtest = ...
                obj.data.Xtest;

            Ytest = ...
                obj.data.Ytest;

            switch obj.costType

                case 'L2'

                    Ypred = ...
                        obj.network.computeYOut(Xtest);

                    num = ...
                        norm(Ypred-Ytest,'fro');

                    den = ...
                        max(norm(Ytest,'fro'),eps);

                    testError = ...
                        num/den;

                case '-loglikelihood'

                    [~,Ypred] = ...
                        max(...
                        obj.network.computeLastH(Xtest),...
                        [],2);

                    [~,Ytarget] = ...
                        max(Ytest,[],2);

                    testError = ...
                        mean(Ypred ~= Ytarget);

                otherwise

                    msg = ...
                        [obj.costType,...
                        ' is not a valid cost function'];

                    error(msg)

            end

        end

    end

    methods (Access = private)

        function init(obj,cParams)

            obj.costType = ...
                cParams.costType;

            obj.designVariable = ...
                cParams.designVariable;

            obj.network = ...
                cParams.network;

            obj.data = ...
                cParams.data;

        end

        function j = computeCost(obj,yOut,Yb)

            [j,~] = ...
                obj.lossFunction(Yb,yOut);

        end

        function dj = computeGradient(obj,yOut,Yb)

            [~,dLF] = ...
                obj.lossFunction(Yb,yOut);

            dj = ...
                obj.network.backprop(Yb,dLF);

        end

        function [J,gc] = lossFunction(obj,y,yOut)

            type = ...
                obj.costType;

            switch type

                case '-loglikelihood'

                    yp = ...
                        min(...
                        max(yOut,1e-11),...
                        1-1e-11);

                    c = ...
                        sum(...
                        (1-y).*(-log(1-yp)) + ...
                        y.*(-log(yp)),...
                        2);

                    J = ...
                        mean(c);

                    gc = ...
                        (yp-y)./...
                        (yp.*(1-yp));

                case 'L2'

                    yp = ...
                        yOut;

                    c = ...
                        (yp-y).^2;

                    J = ...
                        0.5*mean(sum(c,2));

                    gc = ...
                        yp-y;

                otherwise

                    msg = ...
                        [type,...
                        ' is not a valid cost function'];

                    error(msg)

            end

        end

        function [ord,nBatches] = ...
                computeNumberOfBatchesAndOrder(obj)

            nD = ...
                size(obj.data.Xtrain,1);

            batchSize = ...
                obj.computeBatchSize();

            nBatches = ...
                ceil(nD/batchSize);

            if nBatches == 1 || ...
                    nBatches == 0

                ord = ...
                    1:nD;

                nBatches = ...
                    1;

            else

                ord = ...
                    randperm(nD,nD);

            end

            obj.order = ...
                ord;

            obj.nBatches = ...
                nBatches;

        end

        function updateCurrentBatch(obj,moveBatch)

            if moveBatch

                obj.currentBatchIndices = ...
                    obj.computeBatchIndices(...
                    obj.iBatch);

                obj.currentBatchDepleted = ...
                    obj.iBatch == ...
                    obj.nBatches;

                obj.iBatch = ...
                    obj.updateBatchCounter(...
                    obj.iBatch,true);

            elseif isempty(obj.currentBatchIndices)

                obj.currentBatchIndices = ...
                    obj.computeBatchIndices(...
                    obj.iBatch);

                obj.currentBatchDepleted = ...
                    obj.iBatch == ...
                    obj.nBatches;

            end

        end

        function indices = ...
                computeBatchIndices(obj,iBatch)

            batchSize = ...
                obj.computeBatchSize();

            iStart = ...
                (iBatch-1)*batchSize + 1;

            iEnd = ...
                min(...
                iBatch*batchSize,...
                size(obj.data.Xtrain,1));

            indices = ...
                obj.order(iStart:iEnd);

        end

        function iB = ...
                updateBatchCounter(...
                obj,iB,moveBatch)

            if ~moveBatch

                return

            end

            if iB < obj.nBatches

                iB = ...
                    iB + 1;

            elseif iB == obj.nBatches

                iB = ...
                    1;

                nD = ...
                    size(obj.data.Xtrain,1);

                obj.order = ...
                    randperm(nD);

            end

        end

        function batchSize = computeBatchSize(obj)

            if size(obj.data.Xtrain,1) > 200

                batchSize = ...
                    200;

            else

                batchSize = ...
                    size(obj.data.Xtrain,1);

            end

        end

    end

end