classdef SGD < Trainer

    properties (Access = protected)
       learningRate  
    end

    properties (Access = private)
       thetaLowest
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
    end

    methods (Access = public)

        function obj = SGD(s)
            obj.init(s)
            obj.plotter = s.plotter;
            obj.learningRate = s.learningRate;
            obj.maxFunEvals  = 5000;
            obj.optTolerance = 10^(-8);
            obj.timeStop    = Inf([1,1]);
            obj.fvStop      = -Inf;%10^(-4);
            obj.nPlot       = 1;
            obj.MaxEpochs   = 1000;
            obj.earlyStop   = obj.MaxEpochs;
            obj.svepoch     = 0;
            obj.lSearchtype  = 'static';
        end
        
        function train(obj)
           tic
           x0  = obj.designVariable.thetavec;
           F = @(theta,order,i) obj.objectiveFunction.computeCost(theta,order,i); 
           obj.optimize(F,x0);
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
        
        function optimize(obj,F,th0)
            [nD, nB, ~] = obj.objectiveFunction.getBatchSize();
            epsilon0      = obj.learningRate;
            epoch         =  1;iter = -1; funcount =  0; fv = 1;
            alarm         =  0; gnorm = 1; min_testError = 1;
            criteria(1)   = epoch <= obj.MaxEpochs; 
            criteria(2)   = alarm < obj.earlyStop; 
            criteria(3)   = gnorm > obj.optTolerance;
            criteria(4)   = toc < obj.timeStop;
            criteria(5)   = fv > obj.fvStop;
            while all(criteria == 1)

                if nB == 1 || nB == 0
                    order = 1:nD;
                    nB = 1;
                else
                    order = randperm(nD,nD);
                end
                for i = 1:nB
                    if iter == -1
                        th      = th0;
                        epsilon = epsilon0;
                        state   = 'init';   
                    else
                        state   = 'iter';
                    end
                    [f,grad]              = F(th,order, i);  
                    [epsilon,th,funcount] = obj.lineSearch(th,grad,F,f,epsilon,epsilon0,funcount,order,i);                
                    gnorm                 = norm(grad,2);
                    opt.epsilon           = epsilon*gnorm; 
                    opt.gnorm             = gnorm;
                    funcount              = funcount + 1;
                    iter                  = iter + 1;
                    obj.displayIter(epoch,iter,funcount,th,f,opt,state);
                end
                [alarm,min_testError] = obj.validateES(alarm,min_testError);
                epoch = epoch + 1;
                criteria(1)   = epoch <= obj.MaxEpochs; 
                criteria(2)   = alarm < obj.earlyStop; 
                criteria(3)   = gnorm > obj.optTolerance;
                criteria(4)   = toc < obj.timeStop; 
                criteria(5)   = f > obj.fvStop;
            end
            if criteria(1) == 0
                fprintf('Minimization terminated, maximum number of epochs reached %d\n',epoch)
            elseif criteria(2) == 0
                fprintf('Minimization terminated, validation set did not decrease in %d epochs\n',obj.earlyStop)
            elseif criteria(3) == 0
                fprintf('Minimization terminated, reached the optimalty tolerance of %f\n',obj.optTolerance)
            elseif criteria(4) == 0
                fprintf('Minimization terminated, reached the limit time of %f\n',obj.timeStop)
            elseif criteria(5) == 0
                fprintf('Minimization terminated, reached the target function value of %f\n',obj.fvStop)
            else
                fprintf('The operation terminated excepcionally\n')
            end
            %F(th,Xb,Yb);
            %F(th,order, i);
            %th
        end

        function [e,x,funcount] = lineSearch(obj,x,grad,F,fOld,e,e0,funcount,order,i)
            type = obj.lSearchtype;
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
                        [f,~] = F(xnew,order,i);
                        e = e/2;
                        funcount = funcount + 1;
                    end
                    e = 5*e;                
                case 'fminbnd'
                    xnew = @(e1) obj.step(x,e1,grad);
                    f = @(e1) F(xnew(e1),order,i);
                    [e,~] = fminbnd(f,e/10,e*10);
                    xnew = xnew(e);
            end
            x = xnew;
        end 

        function [alarm,min_testError] = validateES(obj,alarm,min_testError)
            [Xtest, Ytest] = obj.objectiveFunction.getTestData();
            [~,y_pred]   = max(obj.objectiveFunction.getOutput(Xtest),[],2);
            [~,y_target] = max(Ytest,[],2);
            testError    = mean(y_pred ~= y_target);
            if testError < min_testError
                obj.thetaLowest = obj.designVariable.thetavec;
                min_testError = testError;
                alarm = 0;
            elseif testError == min_testError
                alarm = alarm + 0.5;
            else
                alarm = alarm + 1;
            end
        end

        function displayIter(obj,epoch,iter,funcount,x,f,opt,state)
            [nD, ~, batchSize] = obj.objectiveFunction.getBatchSize();
            obj.printValues(epoch,funcount,opt,f,iter)
            if obj.isDisplayed == true
                if iter*batchSize==(epoch-1)*nD || iter == -1
                    obj.storeValues(x,f,state,opt);
                    obj.plotMinimization(epoch);
                end
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