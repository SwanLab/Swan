classdef Optimizer_Constrained < Optimizer
    properties
        fhtri
        niter = 0
        optimizer
        maxiter        
    end
    
    methods
        function obj = Optimizer_Constrained(settings,mesh,monitoring)
            obj@Optimizer(settings);
            obj.maxiter = settings.maxiter;
            obj.plotting = settings.plotting;
            obj.printing = settings.printing;
            obj.monitoring = Monitoring(settings,monitoring);
            obj.optimizer = settings.optimizer;
            obj.mesh = mesh;            
        end
        
        function x = solveProblem(obj,x_ini,cost,constraint,istep,nstep)
            cost.computef(x_ini);
            constraint.computef(x_ini);
            obj.plotX(x_ini)
            obj.print(x_ini,obj.niter);
            while obj.stop_criteria && obj.niter < obj.maxiter
                obj.niter = obj.niter+1;
                x = obj.updateX(x_ini,cost,constraint);
                obj.plotX(x)
                obj.print(x,obj.niter);
                obj.monitoring.display(obj.niter,cost,constraint,obj.stop_vars,obj.stop_criteria && obj.niter < obj.maxiter,istep,nstep);
                x_ini = x;
            end
            obj.stop_criteria = 1;
        end
    end
    
%     methods (Abstract)
%         setContraintsCase(obj);
%     end
end
