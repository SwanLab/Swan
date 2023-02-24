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

    methods(Access = public)

        function self = SGD(s)
            self.init(s)
            self.data         = s{1}.data;
            self.learningRate = s{3};
            
            if length(s) <= 4
                self.optTolerance = 10^(-6);
                if size(s{1}.data.Xtrain,1) > 200
                    self.batchSize    = 200;
                else
                    self.batchSize    = size(s{1}.data.Xtrain,1);
                end
                self.MaxEpochs   = 2000*self.batchSize/size(s{1}.data.Xtrain,1);
                self.earlyStop   = self.MaxEpochs;
                self.MaxFunEvals = 2000;
                self.timeStop    = Inf([1,1]);
                self.fvStop      = 10^-4;
                self.nPlot       = 1;
                self.lSearchtype  = 'static';
            else 
                self.optTolerance = s{6}.optTolerance;
                if s{5} > size(s{1}.data.Xtrain,1)
                    self.batchSize = size(s{1}.data.Xtrain,1);
                else
                    self.batchSize = s{5};
                end
                self.MaxEpochs   = s{6}.maxepochs;
                self.earlyStop   = s{6}.earlyStop;
                self.MaxFunEvals = s{6}.maxevals;
                self.timeStop    = s{6}.time;
                self.fvStop      = s{6}.fv;
                self.nPlot       = 1;
                self.lSearchtype  = s{7};
            end
            if length(s) == 8
                self.nPlot = s{8};
            end
        end
        
        function train(self)
           tic
           x0  = self.network.thetavec; 
           F = @(theta,X,Y) self.costFunction(theta,X,Y); 
           self.optimize(F,x0);
           toc
        end 
    end
    
    methods(Access = private)  
        
        function optimize(self,F,th0)
            nD            = size(self.data.Xtrain,1);
            epsilon0      = self.learningRate;
            epoch         =  1;iter = -1; funcount =  0; fv = 1;
            alarm         =  0; gnorm = 1; min_testError = 1;
            nB            = fix(nD/self.batchSize);
            criteria(1)   = epoch <= self.MaxEpochs; 
            criteria(2)   = alarm < self.earlyStop; 
            criteria(3)   = gnorm > self.optTolerance;
            criteria(4)   = toc < self.timeStop;
            criteria(5)   = fv > self.fvStop;
            while all(criteria == 1)
                if nB == 1 || nB == 0
                    order = 1:nD;
                    nB = 1;
                else
                    order = randperm(nD,nD);
                end
                for i = 1:nB
                    [Xb,Yb] = self.createMinibatch(order,i);
                    if iter == -1
                        th      = th0;
                        epsilon = epsilon0;
                        state   = 'init';   
                    else
                        state   = 'iter';
                    end
                    [f,grad]              = F(th,Xb,Yb);  
                    [epsilon,th,funcount] = self.lineSearch(th,grad,F,f,epsilon,epsilon0,funcount,Xb,Yb);                
                    gnorm                 = norm(grad,2);
                    opt.epsilon           = epsilon*gnorm; 
                    opt.gnorm             = gnorm;
                    funcount              = funcount + 1;
                    iter                  = iter + 1;
                    self.displayIter(epoch,iter,funcount,th,f,opt,state);
                end
                [alarm,min_testError] = self.validateES(alarm,min_testError);
                epoch = epoch + 1;
                criteria(1)   = epoch <= self.MaxEpochs; 
                criteria(2)   = alarm < self.earlyStop; 
                criteria(3)   = gnorm > self.optTolerance;
                criteria(4)   = toc < self.timeStop; 
                criteria(5)   = self.network.cost > self.fvStop;
            end
            if criteria(1) == 0
                fprintf('Minimization terminated, maximum number of epochs reached %d\n',epoch)
            elseif criteria(2) == 0
                fprintf('Minimization terminated, validation set did not decrease in %d epochs\n',self.earlyStop)
            elseif criteria(3) == 0
                fprintf('Minimization terminated, reached the optimalty tolerance of %f\n',self.optTolerance)
            elseif criteria(4) == 0
                fprintf('Minimization terminated, reached the limit time of %f\n',self.timeStop)
            elseif criteria(5) == 0
                fprintf('Minimization terminated, reached the target function value of %f\n',self.fvStop)
            else
                fprintf('The operation terminated excepcionally\n')
            end
            F(th,Xb,Yb);
        end

        function [e,x,funcount] = lineSearch(self,x,grad,F,fOld,e,e0,funcount,Xb,Yb)
            type = self.lSearchtype;
            switch type
                case 'static'
                    xnew = self.step(x,e,grad,F,Xb,Yb);
                case 'decay'
                    tau = 50;
                    xnew = self.step(x,e,grad,F,Xb,Yb);
                    e = e - 0.99*e0*30/tau;
                case 'dynamic'
                    f = fOld;
                    xnew = x;
                    while f >= 1.001*(fOld - e*(grad*grad'))
                        xnew = self.step(x,e,grad,F,Xb,Yb);
                        [f,~] = F(xnew,Xb,Yb);
                        e = e/2;
                        funcount = funcount + 1;
                    end
                    e = 5*e;                
                case 'fminbnd'
                    xnew = @(e1) self.step(x,e1,grad,F,Xb,Yb);
                    f = @(e1) F(xnew(e1),Xb,Yb);
                    [e,~] = fminbnd(f,e/10,e*10);
                    xnew = xnew(e);
            end
            x = xnew;
        end 

        function [alarm,min_testError] = validateES(self,alarm,min_testError)
            [~,y_pred]   = max(self.network.getOutput(self.data.Xtest),[],2);
            [~,y_target] = max(self.data.Ytest,[],2);
            testError    = mean(y_pred ~= y_target);
            if testError < min_testError
                self.thetaLowest = self.network.thetavec;
                min_testError = testError;
                alarm = 0;
            elseif testError == min_testError
                alarm = alarm + 0.5;
            else
                alarm = alarm + 1;
            end
        end

        function displayIter(self,epoch,iter,funcount,x,f,opt,state)
            self.printValues(epoch,funcount,opt,f,iter)
            if self.isDisplayed == true
                if iter*self.batchSize==(epoch-1)*length(self.data.Xtrain) || iter == -1
                    self.storeValues(x,f,state,opt);
                    self.plotMinimization(epoch);
                end
            end
        end  

        function printValues(self,epoch,funcount,opt,f,iter)
            formatstr = '%5.0f    %5.0f       %5.0f    %13.6g  %13.6g   %12.3g\n';
            if mod(iter,20) == 0
                fprintf(['                                                        First-order \n', ...
                    'Epoch Iteration  Func-count       f(x)        Step-size       optimality\n']);
            end
            fprintf(formatstr,epoch,iter,funcount,f,opt.epsilon,opt.gnorm);
        end

        function [x,y] = createMinibatch(self,order,i)
               I = self.batchSize;            
               X = self.data.Xtrain;
               Y =  self.data.Ytrain;
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
        function x = step(self,x,e,grad,F,Xb,Yb)
            x = x - e*grad;
        end
    end
end