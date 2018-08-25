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
            obj.printing = settings.printing;
            obj.monitoring = Monitoring.create(settings,mesh,monitoring,settings.plotting);
            obj.optimizer = settings.optimizer;
            obj.mesh = mesh;
        end
        
        function x = solveProblem(obj,x_ini,cost,constraint,istep,nstep)
            x = x_ini;
            cost.computef(x_ini);
            constraint.computef(x_ini);
            obj.monitoring.plotX(x_ini)
            obj.print(x_ini,obj.niter);
            while ~obj.stop_updating && obj.niter < obj.maxiter
                obj.niter = obj.niter+1;
                x = obj.updateX(x_ini,cost,constraint);
                obj.monitoring.refresh(x,obj.niter,cost,constraint,obj.stop_vars,obj.stop_updating && obj.niter < obj.maxiter,istep,nstep);
                obj.print(x,obj.niter);
                obj.writeToFile(istep,cost,constraint)
                x_ini = x;
            end
            obj.printFinal(x);
            % !!????? NEEDED ??????
            obj.stop_updating = 0;
        end
        function printFinal(obj,x)
            if obj.monitoring.plotting_ON
                if obj.printing
                    obj.print(x,obj.niter);
                else
                    obj.printing = 1;
                    obj.print(x,obj.niter);
                    obj.printing = 0;
                end
            end
        end
    end
    
%     methods (Abstract)
%         setContraintsCase(obj);
%     end
end
