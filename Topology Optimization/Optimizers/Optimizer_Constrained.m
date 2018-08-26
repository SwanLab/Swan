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
            x_ini = obj.compute_initial_value(x_ini,cost,constraint);
            cost.computef(x_ini);
            constraint.computef(x_ini);
            obj.plotX(x_ini)
            obj.print(x_ini,obj.niter);
            x = x_ini;
            while ~obj.has_converged && obj.niter < obj.maxiter
                obj.niter = obj.niter+1;
                x = obj.updateX(x_ini,cost,constraint);
                obj.plotX(x)
                obj.print(x,obj.niter);
                obj.writeToFile(istep,cost,constraint)
                obj.monitoring.display(obj.niter,cost,constraint,obj.stop_vars,obj.has_converged && obj.niter < obj.maxiter,istep,nstep);
                x_ini = x;
            end
            obj.has_converged = false;
        end 
        
        function x_ini = compute_initial_value(obj,x_ini,cost,constraint)
            
        end
        
    end
    
%     methods (Abstract)
%         setContraintsCase(obj);
%     end
end
