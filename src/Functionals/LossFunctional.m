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
        function t = getTitleToPlot(~)

            t = 'Data loss';

        end


        function [j,dj] = computeFunctionAndGradient(obj,x)

            %
            % Full-batch value and gradient.
            %

            obj.designVariable.thetavec = x;

            Xb = obj.data.Xtrain;
            Yb = obj.data.Ytrain;

            yOut = obj.network.computeYOut(Xb);

            j  = obj.computeCost(yOut,Yb);
            dj = obj.computeGradient(yOut,Yb);

        end


        function [j,dj] = computeCostAndGradient(obj,x)

            %
            % Compatibility interface expected by CostNN.
            %
            % Keep computeFunctionAndGradient as the actual
            % implementation so existing code remains valid.
            %

            [j,dj] = ...
                obj.computeFunctionAndGradient(x);

        end


        function [j,dj,isBD] = ...
                computeStochasticCostAndGradient(obj,x,moveBatch)

            obj.designVariable.thetavec = x;

            Xt = obj.data.Xtrain;
            Yt = obj.data.Ytrain;

            %
            % Evaluate CURRENT batch.
            %
            [Xb,Yb] = ...
                obj.updateSampledDataSet( ...
                    Xt,Yt,obj.iBatch);

            yOut = ...
                obj.network.computeYOut(Xb);

            j  = obj.computeCost(yOut,Yb);
            dj = obj.computeGradient(yOut,Yb);


            %
            % IMPORTANT:
            %
            % Check whether the batch that has JUST been
            % evaluated is the last batch of the epoch.
            %
            isBD = ...
                obj.isBatchDepleted( ...
                    obj.iBatch,moveBatch);


            %
            % Only after evaluating and checking the current
            % batch do we move the batch counter.
            %
            obj.iBatch = ...
                obj.updateBatchCounter( ...
                    obj.iBatch,moveBatch);

        end


        function testError = getTestError(obj)

            %
            % Regression validation error.
            %
            % The validation metric is the same value-based
            % loss used during training.
            %

            Xtest = obj.data.Xtest;
            Ytest = obj.data.Ytest;

            if isempty(Xtest)

                testError = NaN;
                return

            end

            Ypred = ...
                obj.network.computeYOut(Xtest);

            testError = ...
                obj.computeCost(Ypred,Ytest);

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

            %
            % lossFunction returns the derivative of the
            % per-sample loss with respect to network output.
            %
            % Network.backprop performs the mean over the
            % samples of the current batch.
            %

            [~,dLF] = ...
                obj.lossFunction(Yb,yOut);

            dj = ...
                obj.network.backprop(Yb,dLF);

        end


        function [J,gc] = lossFunction(obj,y,yOut)

            type = obj.costType;

            switch type

                case 'L2'

                    %
                    % Regression loss:
                    %
                    %            1
                    % J = ---------------- sum_i ||yhat_i-y_i||^2
                    %          2 m
                    %
                    % Network.backprop applies the factor 1/m.
                    %

                    e = yOut - y;

                    J = ...
                        0.5 * mean( ...
                            sum(e.^2,2));

                    %
                    % Per-sample derivative before averaging:
                    %
                    % dL/dyhat = yhat - y
                    %
                    gc = e;


                case '-loglikelihood'

                    %
                    % Binary cross-entropy.
                    %
                    % Not used by the constitutive surrogate,
                    % but retained for compatibility.
                    %

                    epsProb = 1e-11;

                    yp = ...
                        min( ...
                            max(yOut,epsProb), ...
                            1-epsProb);

                    c = ...
                        (1-y).*(-log(1-yp)) ...
                        + ...
                        y.*(-log(yp));

                    J = ...
                        mean(sum(c,2));

                    gc = ...
                        (yp-y) ./ ...
                        (yp.*(1-yp));


                otherwise

                    error( ...
                        'LossFunctional:InvalidLoss', ...
                        '%s is not a valid loss function.', ...
                        type);

            end

        end


        function [ord,nBatches] = ...
                computeNumberOfBatchesAndOrder(obj)

            nD = ...
                size(obj.data.Xtrain,1);

            batchSize = ...
                obj.computeBatchSize();

            %
            % We preserve the original batch partition:
            %
            % the remainder is included in the final batch.
            %
            nBatches = ...
                fix(nD/batchSize);

            if nBatches == 1 || nBatches == 0

                ord = 1:nD;
                nBatches = 1;

            else

                %
                % Random permutation for the first epoch.
                %
                ord = randperm(nD);

            end

            obj.order = ord;
            obj.nBatches = nBatches;

        end


        function itIs = ...
                isBatchDepleted(obj,iBatch,moveBatch)

            %
            % True only after evaluating the last batch of
            % the current epoch.
            %

            itIs = ...
                moveBatch ...
                && ...
                iBatch == obj.nBatches;

        end


        function iB = ...
                updateBatchCounter(obj,iB,moveBatch)

            if ~moveBatch
                return
            end


            if iB < obj.nBatches

                %
                % Continue inside the current epoch.
                %
                iB = iB + 1;

            else

                %
                % Current epoch is finished.
                %
                % Start the next epoch at batch 1.
                %
                iB = 1;

                %
                % Random reshuffling:
                % every new epoch receives a new permutation
                % of the training samples.
                %
                obj.reshuffle();

            end

        end


        function reshuffle(obj)

            %
            % With a single full batch, the ordering of the
            % samples has no influence on the batch gradient.
            %

            if obj.nBatches <= 1
                return
            end

            nD = ...
                size(obj.data.Xtrain,1);

            obj.order = ...
                randperm(nD);

        end


        function batchSize = computeBatchSize(obj)

            %
            % Preserve the original batch-size policy.
            %

            if size(obj.data.Xtrain,1) > 200

                batchSize = 200;

            else

                batchSize = ...
                    size(obj.data.Xtrain,1);

            end

        end


        function [x,y] = ...
                updateSampledDataSet( ...
                    obj,Xl,Yl,iBatch)

            %
            % Preserve the original convention in which the
            % remainder is absorbed by the final batch.
            %

            iB = iBatch;

            batchSize = ...
                obj.computeBatchSize();

            nB = batchSize;


            if iB == fix(size(Xl,1)/nB)

                plus = ...
                    mod(size(Xl,1),nB);

                x = ...
                    zeros( ...
                        [nB+plus,size(Xl,2)]);

                y = ...
                    zeros( ...
                        [nB+plus,size(Yl,2)]);

            else

                plus = 0;

                x = ...
                    zeros( ...
                        [nB,size(Xl,2)]);

                y = ...
                    zeros( ...
                        [nB,size(Yl,2)]);

            end


            cont = 1;

            for jB = ...
                    (iB-1)*nB+1 : ...
                    iB*nB+plus

                x(cont,:) = ...
                    Xl(obj.order(jB),:);

                y(cont,:) = ...
                    Yl(obj.order(jB),:);

                cont = cont + 1;

            end

        end

    end

end