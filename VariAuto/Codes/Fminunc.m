classdef Fminunc < Trainer

    properties (Access = private)
        opt
        
    end

    methods(Access = public)
        function obj = Fminunc(s)
            obj.init(s);
            s.optTolerance = 1*10^-10;
            s.maxevals = 5000;
            s.nPlot = 1;
            obj.opt   = obj.setSolverOptions(s);
            obj.nPlot = s.nPlot;
            obj.Xtrain  = s.Xtrain;
            obj.Ytrain  = s.Ytrain;
        end

        function train(obj)
            x0  = obj.designVariable.thetavec;
            c   = obj.costFunction;
            F = @(theta) obj.costFunction.computeCost(theta,obj.Xtrain,obj.Ytrain);
            fminunc(F,x0,obj.opt); 
        end
    end

    methods(Access = private)

        function opt = setSolverOptions(obj,s)
           opt = optimoptions(@fminunc);
           opt.SpecifyObjectiveGradient = false;
           opt.Algorithm                = 'quasi-newton';
           opt.OptimalityTolerance      = s.optTolerance;
           opt.MaxIterations            = s.maxevals*5;
           opt.MaxFunctionEvaluations   = s.maxevals; 
           if obj.isDisplayed == true
                args = [];
                opt.Display        = 'iter';
                opt.CheckGradients = true;
                opt.OutputFcn      = @(theta,optimvalues,state)obj.myoutput(theta,optimvalues,state,args);
           end
        end 

        function stop = myoutput(obj,x,optimvalues,state,args)
            stop         = false;
            f            = optimvalues.fval;
            opti.epsilon = optimvalues.stepsize;
            opti.gnorm   = optimvalues.firstorderopt;
            iter         = optimvalues.iteration;
            if iter == 0
                opti.epsilon = 1;
            end
            obj.storeValues(x,f,state,opti);
            obj.plotMinimization(iter);                                
        end
    end
end