classdef SGD < Trainer

    properties (Access = protected)
       lSearchtype
       learningRate
       MaxEpochs
       MaxFunEvals
       optTolerance
       earlyStop
       timeStop
       fvStop
       batchSize
       thetaLowest
    end

    methods (Access = public)

        function obj = SGD(s)
            obj.init(s)
            obj.learningRate = s{3};
            
            if length(s) <= 4
                obj.optTolerance = 10^(-6);
                if size(obj.CostFunction.data.Xtrain,1) > 200
                    obj.batchSize    = 200;
                else
                    obj.batchSize    = size(s{1}.data.Xtrain,1);
                end
                obj.MaxEpochs   = 2000*obj.batchSize/size(obj.CostFunction.data.Xtrain,1);
                obj.earlyStop   = obj.MaxEpochs;
                obj.MaxFunEvals = 2000;
                obj.timeStop    = Inf([1,1]);
                obj.fvStop      = 10^-4;
                obj.nPlot       = 1;
                obj.lSearchtype  = 'static';
            else 
                obj.optTolerance = s{6}.optTolerance;
                if s{5} > size(obj.CostFunction.data.Xtrain,1)
                    obj.batchSize = size(obj.CostFunction.data.Xtrain,1);
                else
                    obj.batchSize = s{5};
                end
                obj.MaxEpochs   = s{6}.maxepochs;
                obj.earlyStop   = s{6}.earlyStop;
                obj.MaxFunEvals = s{6}.maxevals;
                obj.timeStop    = s{6}.time;
                obj.fvStop      = s{6}.fv;
                obj.nPlot       = 1;
                obj.lSearchtype  = s{7};
            end
            if length(s) == 8
                obj.nPlot = s{8};
            end
        end
        
        function train(obj)
           tic
           x0  = obj.CostFunction.designVariable.thetavec; 
           F = @(theta,X,Y) obj.costFunction(theta,X,Y); 
           obj.optimize(F,x0);
           toc
        end 
    end
    
    methods(Access = private)  
        
        function optimize(obj,F,th0)
            nD            = size(obj.CostFunction.data.Xtrain,1);
            epsilon0      = obj.learningRate;
            epoch         =  1;iter = -1; funcount =  0; fv = 1;
            alarm         =  0; gnorm = 1; min_testError = 1;
            nB            = fix(nD/obj.batchSize);
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
                    [Xb,Yb] = obj.createMinibatch(order,i);
                    if iter == -1
                        th      = th0;
                        epsilon = epsilon0;
                        state   = 'init';   
                    else
                        state   = 'iter';
                    end
                    [f,grad]              = F(th,Xb,Yb);  
                    [epsilon,th,funcount] = obj.lineSearch(th,grad,F,f,epsilon,epsilon0,funcount,Xb,Yb);                
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
                criteria(5)   = obj.CostFunction.cost > obj.fvStop;
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
            F(th,Xb,Yb);
        end

        function [e,x,funcount] = lineSearch(obj,x,grad,F,fOld,e,e0,funcount,Xb,Yb)
            type = obj.lSearchtype;
            switch type
                case 'static'
                    xnew = obj.step(x,e,grad,F,Xb,Yb);
                case 'decay'
                    tau = 50;
                    xnew = obj.step(x,e,grad,F,Xb,Yb);
                    e = e - 0.99*e0*30/tau;
                case 'dynamic'
                    f = fOld;
                    xnew = x;
                    while f >= 1.001*(fOld - e*(grad*grad'))
                        xnew = obj.step(x,e,grad,F,Xb,Yb);
                        [f,~] = F(xnew,Xb,Yb);
                        e = e/2;
                        funcount = funcount + 1;
                    end
                    e = 5*e;                
                case 'fminbnd'
                    xnew = @(e1) obj.step(x,e1,grad,F,Xb,Yb);
                    f = @(e1) F(xnew(e1),Xb,Yb);
                    [e,~] = fminbnd(f,e/10,e*10);
                    xnew = xnew(e);
            end
            x = xnew;
        end 

        function [alarm,min_testError] = validateES(obj,alarm,min_testError)
            [~,y_pred]   = max(obj.optimizationProblem.getOutput(obj.data.Xtest),[],2);
            [~,y_target] = max(obj.data.Ytest,[],2);
            testError    = mean(y_pred ~= y_target);
            if testError < min_testError
                obj.thetaLowest = obj.optimizationProblem.designVariable.thetavec;
                min_testError = testError;
                alarm = 0;
            elseif testError == min_testError
                alarm = alarm + 0.5;
            else
                alarm = alarm + 1;
            end
        end

        function displayIter(obj,epoch,iter,funcount,x,f,opt,state)
            obj.printValues(epoch,funcount,opt,f,iter)
            if obj.isDisplayed == true
                if iter*obj.batchSize==(epoch-1)*length(obj.CostFunction.data.Xtrain) || iter == -1
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
        end

        function [x,y] = createMinibatch(obj,order,i)
               I = obj.batchSize;            
               X = obj.CostFunction.data.Xtrain;
               Y =  obj.data.Ytrain;
               cont = 1;
               if i == fix(size(X,1)/I)
                   plus = mod(size(X,1),I);
                   x = zeros([I+plus,size(X,2)]);
                   y = zeros([I+plus,size(Y,2)]);
               else
                   plus = 0;
                   x = zeros([I,size(X,2)]);
                   y = zeros([I,size(Y,2)]);
               end
               for j = (i-1)*I+1:i*I+plus
                   x(cont,:) = X(order(j),:);
                   y(cont,:) = Y(order(j),:);
                   cont = cont+1;
               end
           end
    end
    methods (Access = protected)
        function x = step(obj,x,e,grad,F,Xb,Yb)
            x = x - e*grad;
        end
    end
end