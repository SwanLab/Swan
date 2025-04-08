classdef SGD < Trainer

    properties (Access = private)
       fvStop
       lSearchType
       MaxEpochs
       maxFunEvals
       optTolerance
       earlyStop
       timeStop
       plotter
       svepoch
       fplot
       learningRate  

    end

    methods (Access = public)

        function obj = SGD(s)
            obj.init(s)
            obj.plotter = s.plotter;
            obj.learningRate = s.learningRate;
            obj.maxFunEvals  = 5000;
            obj.optTolerance = 1e-8;
            obj.timeStop    = Inf([1,1]);
            obj.fvStop      = 1e-4;
            obj.nPlot       = 1;
            obj.MaxEpochs   = 1000;
            obj.earlyStop   = obj.MaxEpochs;
            obj.svepoch     = 0;
            obj.lSearchType = 'static';
        end
        
        function compute(obj)
           tic
           x0  = obj.designVariable.thetavec;
           obj.optimize(x0);
           toc
        end

        function plotCostFunc(obj)
            figure(3);
            % epoch = 1:obj.MaxEpochs;
            epoch = 1:length(obj.fplot);
            % plot(epoch,obj.fplot,'-o');
            grid on
            loglog(epoch,obj.fplot,'LineWidth',1.8);
            xlabel('Epochs')
            ylabel('Function Values')
            title('Cost Function')
            xlim([1,inf])
        end

    end
    
    methods(Access = private)  
        
        function optimize(obj,th0)
            epsilon      = obj.learningRate;
            iter         = -1;
            funcount     = 0;
            alarm        = 0;
            minTestError = 1;
            KPI.epoch = 1;
            KPI.alarm = 0;
            KPI.gnorm = 1;
            KPI.cost  = 1;

            while obj.isCriteriaMet(KPI) == false
                state   = 'iter';
                if iter == -1
                    theta   = th0;
                    state   = 'init';   
                end
                newEpoch  = true;
                moveBatch = true;
                while obj.objectiveFunction.isBatchDepleted == false || newEpoch
                    [f, grad] = obj.computeStochasticFunctionAndGradient(theta, moveBatch);
                    [epsilon,theta,funcount] = obj.lineSearch(theta,grad,f,epsilon,funcount);

                    funcount  = funcount + 1;
                    iter      = iter + 1;
                    newEpoch  = false;
                    
                    KPI.cost  = f;
                    KPI.gnorm = norm(grad,2);
                    obj.displayIter(iter,funcount,theta,epsilon,state,KPI);
                end
                KPI.epoch = KPI.epoch + 1;
                [KPI.alarm,minTestError] = obj.objectiveFunction.validateES(alarm,minTestError);
            end
        end

        function [f, grad] = computeStochasticFunctionAndGradient(obj, theta, moveBatch)
            obj.objectiveFunction.setBatchMover(moveBatch);
            obj.objectiveFunction.computeStochasticFunctionAndGradient(theta);
            f    = obj.objectiveFunction.value;
            grad = obj.objectiveFunction.gradient;
        end

        function [e,x,funcount] = lineSearch(obj,x,grad,fOld,e,funcount)
            moveBatch = false;
            F = @(theta) obj.computeStochasticFunctionAndGradient(theta, moveBatch);
            type = obj.lSearchType;
            switch type
                case 'static'
                    xnew = obj.step(x,e,grad);
                case 'decay'
                    xnew = obj.step(x,e,grad);
                    e = e * (1 - 1e-3);
                case 'dynamic'
                    f = fOld;
                    xnew = x;
                    while f >= 1.001*(fOld - e*(grad*grad'))
                        xnew = obj.step(x,e,grad);
                        [f, ~] = F(xnew);
                        e = e/2;
                        funcount = funcount + 1;
                    end
                    e = 5*e;                
                case 'fminbnd'
                    xnew = @(e1) obj.step(x,e1,grad);
                    f = @(e1) F(xnew(e1));
                    [e,~] = fminbnd(f,e/10,e*10);
                    xnew = xnew(e);
            end
            x = xnew;
        end

        function criteria = updateCriteria(obj, KPI)
            criteria(1)   = KPI.epoch <= obj.MaxEpochs; 
            criteria(2)   = KPI.alarm < obj.earlyStop;
            criteria(3)   = KPI.gnorm > obj.optTolerance;
            criteria(4)   = toc < obj.timeStop;
            criteria(5)   = KPI.cost > obj.fvStop;
        end

        function itIs = isCriteriaMet(obj, KPI)
            criteria = obj.updateCriteria(KPI);
            failedIdx = find(~criteria, 1);
            itIs = false;
            if ~isempty(failedIdx)
                itIs = true;
                msg = {sprintf('Minimization terminated, maximum number of epochs reached %d\n',KPI.epoch), ...
                       sprintf('Minimization terminated, validation set did not decrease in %d epochs\n',obj.earlyStop), ...
                       sprintf('Minimization terminated, reached the optimalty tolerance of %f\n',obj.optTolerance), ...
                       sprintf('Minimization terminated, reached the limit time of %f\n',obj.timeStop), ...
                       sprintf('Minimization terminated, reached the target function value of %f\n',obj.fvStop)};
                if failedIdx >= 1 && failedIdx <= length(msg)
                    fprintf('%s', msg{failedIdx})
                end
            end
        end

        function displayIter(obj,iter,funcount,x,epsilon,state,KPI)
            opt.epsilon = epsilon*KPI.gnorm; 
            opt.gnorm   = KPI.gnorm;
            obj.printValues(KPI.epoch,funcount,opt,KPI.cost,iter)
            if obj.isDisplayed && (~mod(epoch, 25) || iter == -1)
                obj.storeValues(x,f,state,opt);
            end
        end  

        function printValues(obj,epoch,funcount,opt,f,iter)
            formatstr = '%5.0f    %5.0f       %5.0f    %13.6g  %13.6g   %12.3g\n';
            if mod(iter,20) == 0
                fprintf(['                                                        First-order \n', ...
                    'Epoch Iteration  Func-count       f(x)        Step-size       optimality\n']);
            end
            fprintf(formatstr,epoch,iter,funcount,f,opt.epsilon,opt.gnorm);
            if epoch ~= obj.svepoch
                obj.svepoch = epoch;
                obj.fplot(1,epoch) = f;
            end
        end

    end
    
    methods (Access = protected)
        function x = step(obj,x,e,grad)
            x = x - e*grad;
        end
    end
end