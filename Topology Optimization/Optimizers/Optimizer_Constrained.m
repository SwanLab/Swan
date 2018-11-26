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
            obj.print(x_ini,obj.niter);
            x_ini = obj.compute_initial_value(x_ini,cost,constraint); % !! REMOVE WHEN DesginVariableInitializer CONSIDERS Projected_Slerp INITIAL GUESS !!
            x = x_ini;
            cost.computeCostAndGradient(x_ini);
            constraint.computeCostAndGradient(x_ini);
%             obj.monitoring.plotX(x_ini);
            obj.monitoring.refresh(x,obj.niter,cost,constraint,obj.stop_vars,obj.has_converged || obj.niter > obj.maxiter*(istep/nstep),istep,nstep);
            
            while ~obj.has_converged && obj.niter < obj.maxiter*(istep/nstep)
                obj.niter = obj.niter+1;
                x = obj.updateX(x_ini,cost,constraint);
                obj.monitoring.refresh(x,obj.niter,cost,constraint,obj.stop_vars,obj.has_converged || obj.niter > obj.maxiter*(istep/nstep),istep,nstep);
                obj.print(x,obj.niter);
                obj.writeToFile(istep,cost,constraint)
                x_ini = x;
            end
            obj.printFinal(x);
            % !!????? NEEDED ??????
            obj.has_converged = 0;
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
        
        function x_ini = compute_initial_value(obj,x_ini,cost,constraint)
            % !! REMOVE WHEN DesginVariableInitializer CONSIDERS
            % Projected_Slerp INITIAL GUESS !!
        end
        
        
    end
    
%     methods (Abstract)
%         setContraintsCase(obj);
%     end
end
