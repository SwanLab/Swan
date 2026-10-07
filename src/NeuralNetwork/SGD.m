classdef SGD < Trainer

    properties (Access = private)

        fvStop
        lSearchType

        MaxEpochs
        optTolerance
        earlyStop
        timeStop

        plotter
        fplot

        learningRate

        %
        % Timer identifier.
        %
        tStart

        %
        % Best model according to validation loss.
        %
        bestTheta
        bestEpoch
        bestValidationError

    end


    methods (Access = public)

        function obj = SGD(s)

            obj.init(s);

            obj.plotter = ...
                s.plotter;

            obj.learningRate = ...
                s.learningRate;

            obj.MaxEpochs = ...
                s.maxEpochs;


            %
            % Optimization tolerances.
            %
            obj.optTolerance = ...
                1e-8;

            obj.timeStop = ...
                Inf;


            %
            % IMPORTANT:
            %
            % This value now refers to
            %
            % J = 1/(2m) sum_i ||yhat_i-y_i||^2
            %
            % and is therefore not numerically equivalent
            % to the old sqrt(sum(e.^2)) criterion.
            %
            % We preserve the previous value for now so that
            % hyperparameter tuning is not mixed with the
            % structural corrections.
            %
            if isfield(s,'fvStop')
                obj.fvStop = s.fvStop;
            else
                obj.fvStop = -Inf;
            end


            %
            % Validation patience.
            %
            % Preserve old behavior unless explicitly supplied.
            %
            if isfield(s,'earlyStop')

                obj.earlyStop = ...
                    s.earlyStop;

            else

                obj.earlyStop = ...
                    obj.MaxEpochs;

            end


            %
            % Current surrogate training uses a fixed
            % learning rate.
            %
            if isfield(s,'lSearchType')

                obj.lSearchType = ...
                    s.lSearchType;

            else

                obj.lSearchType = ...
                    'static';

            end


            obj.nPlot = 1;

            obj.fplot = [];

            obj.bestTheta = [];
            obj.bestEpoch = 0;
            obj.bestValidationError = Inf;

        end


        function compute(obj)

            %
            % Use an explicit timer identifier so that
            % other tic calls cannot reset this timer.
            %
            obj.tStart = tic;

            x0 = ...
                obj.designVariable.thetavec;

            obj.optimize(x0);

            elapsedTime = ...
                toc(obj.tStart);

            fprintf( ...
                'Training time: %.6f s\n', ...
                elapsedTime);

        end


        function f = getHistory(obj)

            f = ...
                obj.fplot;

        end


        function epoch = getBestEpoch(obj)

            epoch = ...
                obj.bestEpoch;

        end


        function value = getBestValidationError(obj)

            value = ...
                obj.bestValidationError;

        end


        function plotCostFunc(obj)

            figure(3);

            epoch = ...
                1:length(obj.fplot);

            grid on

            loglog( ...
                epoch, ...
                obj.fplot, ...
                'LineWidth',1.8);

            xlabel('Epochs')
            ylabel('Function Values')
            title('Cost Function')

            xlim([1,inf])

        end

    end


    methods (Access = private)

        function optimize(obj,th0)

            %
            % Initial learning rate.
            %
            epsilon = ...
                obj.learningRate;


            %
            % Number of parameter updates.
            %
            iter = 0;


            %
            % Number of objective evaluations.
            %
            funcount = 0;


            %
            % Validation loss is not bounded by one.
            %
            minTestError = ...
                Inf;


            %
            % KPI.epoch means COMPLETED epochs.
            %
            KPI.epoch = 0;

            KPI.alarm = 0;

            %
            % Before the first epoch there is no convergence
            % information.
            %
            KPI.gnorm = ...
                Inf;

            KPI.cost = ...
                Inf;


            theta = ...
                th0;


            %
            % Initial best model.
            %
            % It will only be restored if a finite validation
            % error is actually observed.
            %
            obj.bestTheta = ...
                theta;

            obj.bestEpoch = ...
                0;

            obj.bestValidationError = ...
                Inf;


            while ~obj.isCriteriaMet(KPI)

                % ==================================================
                % MINI-BATCH UPDATES: ONE COMPLETE EPOCH
                % ==================================================

                epochDepleted = ...
                    false;


                while ~epochDepleted

                    %
                    % LossFunctional evaluates the current
                    % mini-batch and then advances its internal
                    % batch counter.
                    %
                    moveBatch = true;

                    [fBatch,gradBatch] = ...
                        obj.computeStochasticFunctionAndGradient( ...
                            theta,moveBatch);


                    %
                    % CostNN reports whether the batch that has
                    % JUST been evaluated was the final one.
                    %
                    epochDepleted = ...
                        obj.objectiveFunction.isBatchDepleted;


                    %
                    % Stop immediately if training diverges
                    % numerically.
                    %
                    if ...
                            ~isfinite(fBatch) ...
                            || ...
                            any(~isfinite(gradBatch(:)))

                        error( ...
                            'SGD:NonFiniteBatch', ...
                            ['Non-finite objective or gradient ', ...
                             'detected during mini-batch training.']);

                    end


                    %
                    % Parameter update.
                    %
                    [epsilon,theta] = ...
                        obj.lineSearch( ...
                            theta, ...
                            gradBatch, ...
                            fBatch, ...
                            epsilon);


                    iter = ...
                        iter + 1;

                    funcount = ...
                        funcount + 1;

                end


                % ==================================================
                % END OF EPOCH
                % ==================================================
                %
                % Convergence criteria must be evaluated using
                % the COMPLETE training objective, not the final
                % mini-batch.
                %

                [fFull,gradFull] = ...
                    obj.computeFunctionAndGradient(theta);

                funcount = ...
                    funcount + 1;


                if ...
                        ~isfinite(fFull) ...
                        || ...
                        any(~isfinite(gradFull(:)))

                    error( ...
                        'SGD:NonFiniteFullBatch', ...
                        ['Non-finite full-batch objective or ', ...
                         'gradient detected after an epoch.']);

                end


                %
                % One complete pass through the training set
                % has been completed.
                %
                KPI.epoch = ...
                    KPI.epoch + 1;


                %
                % First-order information from the complete
                % training problem.
                %
                KPI.cost = ...
                    fFull;

                KPI.gnorm = ...
                    norm(gradFull(:),2);


                % ==================================================
                % VALIDATION / EARLY STOPPING
                % ==================================================

                previousMinTestError = ...
                    minTestError;


                [KPI.alarm,minTestError] = ...
                    obj.objectiveFunction.validateES( ...
                        KPI.alarm, ...
                        minTestError);


                %
                % A strict reduction of minTestError means
                % the current theta is the best validation
                % model observed so far.
                %
                if minTestError < previousMinTestError

                    obj.bestTheta = ...
                        theta;

                    obj.bestEpoch = ...
                        KPI.epoch;

                    obj.bestValidationError = ...
                        minTestError;

                end


                % ==================================================
                % MONITORING
                % ==================================================

                obj.displayIter( ...
                    iter, ...
                    funcount, ...
                    epsilon, ...
                    KPI);

            end


            % ======================================================
            % MODEL SELECTION
            % ======================================================
            %
            % If a validation set was available, early stopping
            % selects the parameters corresponding to the lowest
            % validation loss, not the final epoch.
            %

            if isfinite(obj.bestValidationError)

                theta = ...
                    obj.bestTheta;

                fprintf( ...
                    ['Restoring parameters from epoch %d ', ...
                     '(best validation loss = %.6e).\n'], ...
                    obj.bestEpoch, ...
                    obj.bestValidationError);

            end


            %
            % With no validation set, bestValidationError
            % remains Inf and theta remains the final iterate.
            %
            obj.designVariable.thetavec = ...
                theta;

        end


        function [f,grad] = ...
                computeStochasticFunctionAndGradient( ...
                    obj,theta,moveBatch)

            obj.objectiveFunction. ...
                setBatchMover(moveBatch);

            obj.objectiveFunction. ...
                computeStochasticFunctionAndGradient(theta);

            f = ...
                obj.objectiveFunction.value;

            grad = ...
                obj.objectiveFunction.gradient;

        end


        function [f,grad] = ...
                computeFunctionAndGradient(obj,theta)

            %
            % Complete training-set objective and gradient.
            %

            obj.objectiveFunction. ...
                computeFunctionAndGradient(theta);

            f = ...
                obj.objectiveFunction.value;

            grad = ...
                obj.objectiveFunction.gradient;

        end


        function [e,x] = ...
                lineSearch( ...
                    obj,x,grad,~,e)

            type = ...
                obj.lSearchType;


            switch type

                case 'static'

                    %
                    % Fixed-step mini-batch gradient descent:
                    %
                    % theta_{k+1}
                    % =
                    % theta_k - eta grad J_B(theta_k)
                    %
                    xnew = ...
                        obj.step(x,e,grad);


                case 'decay'

                    %
                    % Preserve the legacy decay rule.
                    %
                    % This is NOT the current training mode.
                    %
                    xnew = ...
                        obj.step(x,e,grad);

                    e = ...
                        e * (1 - 1e-3);


                case {'dynamic','fminbnd'}

                    %
                    % IMPORTANT:
                    %
                    % The current LossFunctional advances its
                    % internal batch immediately after the
                    % gradient evaluation.
                    %
                    % Therefore an additional objective
                    % evaluation here would use another batch,
                    % making Armijo/fminbnd mathematically
                    % inconsistent with the gradient direction.
                    %
                    % These modes must only be enabled after
                    % LossFunctional is modified to cache and
                    % reuse the last evaluated batch.
                    %
                    error( ...
                        'SGD:LineSearchBatchConsistency', ...
                        ['Line-search type "%s" requires ', ...
                         're-evaluation of the same mini-batch. ', ...
                         'Update LossFunctional to cache the ', ...
                         'current batch before enabling this ', ...
                         'line-search mode.'], ...
                        type);


                otherwise

                    error( ...
                        'SGD:LineSearchType', ...
                        'Unknown line-search type "%s".', ...
                        type);

            end


            x = ...
                xnew;

        end


        function criteria = ...
                updateCriteria(obj,KPI)

            %
            % Maximum number of COMPLETED epochs.
            %
            criteria(1) = ...
                KPI.epoch < obj.MaxEpochs;


            %
            % Validation patience.
            %
            criteria(2) = ...
                KPI.alarm < obj.earlyStop;


            %
            % First-order stationarity of the COMPLETE
            % training objective.
            %
            criteria(3) = ...
                KPI.gnorm > obj.optTolerance;


            %
            % Wall-clock limit measured from SGD.compute().
            %
            criteria(4) = ...
                toc(obj.tStart) < obj.timeStop;


            %
            % COMPLETE training objective threshold.
            %
            criteria(5) = ...
                KPI.cost > obj.fvStop;

        end


        function itIs = ...
                isCriteriaMet(obj,KPI)

            criteria = ...
                obj.updateCriteria(KPI);


            failedIdx = ...
                find(~criteria,1);


            itIs = ...
                ~isempty(failedIdx);


            if ~itIs
                return
            end


            msg = { ...
                sprintf( ...
                    ['Minimization terminated: maximum ', ...
                     'number of epochs reached (%d).\n'], ...
                    KPI.epoch), ...
                sprintf( ...
                    ['Minimization terminated: validation ', ...
                     'loss did not improve within the ', ...
                     'allowed patience (%g).\n'], ...
                    obj.earlyStop), ...
                sprintf( ...
                    ['Minimization terminated: full-batch ', ...
                     'gradient norm reached the optimality ', ...
                     'tolerance %.3e.\n'], ...
                    obj.optTolerance), ...
                sprintf( ...
                    ['Minimization terminated: time limit ', ...
                     'reached (%.3g s).\n'], ...
                    obj.timeStop), ...
                sprintf( ...
                    ['Minimization terminated: full-batch ', ...
                     'objective reached the target value ', ...
                     '%.3e.\n'], ...
                    obj.fvStop) ...
                };


            fprintf( ...
                '%s', ...
                msg{failedIdx});

        end


        function displayIter( ...
                obj,iter,funcount,epsilon,KPI)

            %
            % Do NOT label epsilon*||grad|| as the actual
            % mini-batch step size. The optimization uses
            % stochastic gradients, whereas KPI.gnorm is the
            % full-batch gradient norm.
            %
            % Report the learning rate itself instead.
            %

            obj.printValues( ...
                KPI.epoch, ...
                iter, ...
                funcount, ...
                KPI.cost, ...
                epsilon, ...
                KPI.gnorm);

        end


        function printValues( ...
                obj,epoch,iter,funcount,f,epsilon,gnorm)

            formatstr = ...
                ['%5.0f    %9.0f    %10.0f    ', ...
                 '%13.6g    %13.6g    %13.6g\n'];


            if ...
                    epoch == 1 ...
                    || ...
                    ~mod(epoch,20)

                fprintf( ...
                    ['                                                         ', ...
                     'First-order\n', ...
                     'Epoch   Iteration   Func-count        f(x)       ', ...
                     'Learning-rate     optimality\n']);

            end


            fprintf( ...
                formatstr, ...
                epoch, ...
                iter, ...
                funcount, ...
                f, ...
                epsilon, ...
                gnorm);


            %
            % Exactly one complete-training loss per epoch.
            %
            obj.fplot(1,epoch) = ...
                f;

        end

    end


    methods (Access = protected)

        function x = step(obj,x,e,grad)

            %
            % Gradient-descent step.
            %
            x = ...
                x - e*grad;

        end

    end

end