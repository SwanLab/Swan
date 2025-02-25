classdef SGD < Trainer

    properties (Access = private)
       fvStop
       lSearchtype
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
            obj.optTolerance = 10^(-8);
            obj.timeStop    = Inf([1,1]);
            obj.fvStop      = 10^(-4);
            obj.nPlot       = 1;
            obj.MaxEpochs   = 1000;
            obj.earlyStop   = obj.MaxEpochs;
            obj.svepoch     = 0;
            obj.lSearchtype  = 'static';
        end
        
        function compute(obj)
           tic
           x0  = obj.designVariable.thetavec;
           isStochastic = true;
           obj.optimize(x0);
           toc
        end

        function plotCostFunc(obj)
            figure(3);
            % epoch = 1:obj.MaxEpochs;
            epoch = 1:length(obj.fplot);
            % plot(epoch,obj.fplot,'-o');
            plot(epoch,obj.fplot,'LineWidth',1.8);
            xlabel('Epochs')
            ylabel('Function Values')
            title('Cost Function')
            xlim([1,inf])
        end

    end
    
    methods(Access = private)  
        
        function optimize(obj,th0)
            epsilon0      = obj.learningRate;
            epoch         = 1;
            iter          = -1;
            funcount      = 0;
            fv            = 1;
            alarm         = 0;
            gnorm         = 1;
            minTestError = 1;

            

            criteria = obj.updateCriteria(epoch, alarm, gnorm, fv);
            while all(criteria == 1)

                state   = 'iter';
                if iter == -1
                    theta      = th0;
                    epsilon = epsilon0;
                    state   = 'init';   
                end

                batchesDepleted = false;
                moveBatch = true;
                isStochastic = true;
                while batchesDepleted == false
                    obj.objectiveFunction.computeFunctionAndGradient(theta,isStochastic,moveBatch)
                    f = obj.objectiveFunction.value;
                    grad = obj.objectiveFunction.gradient;                    
                    batchesDepleted = obj.objectiveFunction.isBatchDepleted;
                    [epsilon,theta,funcount] = obj.lineSearch(theta,grad,f,epsilon,epsilon0,funcount);                
                    gnorm    = norm(grad,2);
                    funcount = funcount + 1;
                    iter     = iter + 1;
                    obj.displayIter(epoch,iter,funcount,theta,f,gnorm,epsilon,state);
                end
                
                [alarm,minTestError] = obj.objectiveFunction.validateES(alarm,minTestError);
                epoch = epoch + 1;
                criteria = obj.updateCriteria(epoch, alarm, gnorm, f);
            end
        end

        function [e,x,funcount] = lineSearch(obj,x,grad,fOld,e,e0,funcount)
            isStochastic = true;
            F = @(theta,moveBatch) obj.objectiveFunction.computeFunctionAndGradient(theta,isStochastic,moveBatch);

            type = obj.lSearchtype;
            moveBatch = false;
            switch type
                case 'static'
                    xnew = obj.step(x,e,grad);
                case 'decay'
                    tau = 50;
                    xnew = obj.step(x,e,grad);
                    e = e - 0.99*e0*30/tau;
                case 'dynamic'
                    f = fOld;
                    xnew = x;
                    while f >= 1.001*(fOld - e*(grad*grad'))
                        xnew = obj.step(x,e,grad);
                        [f,~,~] = F(xnew, moveBatch);
                        e = e/2;
                        funcount = funcount + 1;
                    end
                    e = 5*e;                
                case 'fminbnd'
                    xnew = @(e1) obj.step(x,e1,grad);
                    f = @(e1) F(xnew(e1), moveBatch);
                    [e,~] = fminbnd(f,e/10,e*10);
                    xnew = xnew(e);
            end
            x = xnew;
        end

        function criteria = updateCriteria(obj, epoch, alarm, gnorm, cost)
            criteria(1)   = epoch <= obj.MaxEpochs; 
            criteria(2)   = alarm < obj.earlyStop; 
            criteria(3)   = gnorm > obj.optTolerance;
            criteria(4)   = toc < obj.timeStop;
            criteria(5)   = cost > obj.fvStop;

            failedIdx = find(~criteria, 1);
            if ~isempty(failedIdx)
                msg = {sprintf('Minimization terminated, maximum number of epochs reached %d\n',epoch), ...
                       sprintf('Minimization terminated, validation set did not decrease in %d epochs\n',obj.earlyStop), ...
                       sprintf('Minimization terminated, reached the optimalty tolerance of %f\n',obj.optTolerance), ...
                       sprintf('Minimization terminated, reached the limit time of %f\n',obj.timeStop), ...
                       sprintf('Minimization terminated, reached the target function value of %f\n',obj.fvStop)};

                if failedIdx >= 1 && failedIdx <= length(msg)
                    fprintf('%s', msg{failedIdx})
                end
            end
        end

        function displayIter(obj,epoch,iter,funcount,x,f,gnorm,epsilon,state)
            opt.epsilon        = epsilon*gnorm; 
            opt.gnorm          = gnorm;
            obj.printValues(epoch,funcount,opt,f,iter)

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