classdef Optimizer_Constrained < Optimizer
    
    properties (Access = public)
        fhtri
        niter = 0
        optimizer
        maxiter
    end
    
    properties (Access = protected)
        monitoring
    end
    
    properties (Access = private)
        printer
        printing
        fileName
        mesh
        printMode
    end
    
    methods (Access = public)
        
        function obj = Optimizer_Constrained(settings,mesh,monitoring)
            obj@Optimizer(settings);
            obj.init(settings,mesh);
            obj.monitoring = Monitoring.create(settings,mesh,monitoring,settings.plotting);
        end
        
        function x = solveProblem(obj,x_ini,cost,constraint,istep,nstep)
            obj.createOptimizerPrinter(cost,constraint);

%            obj.print(x_ini,obj.niter,cost,constraint);
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
                obj.print(x,obj.niter,cost,constraint);
                obj.writeToFile(istep,cost,constraint)
                x_ini = x;
            end
            obj.printFinal(x,cost,constraint);
            % !!????? NEEDED ??????
            obj.has_converged = 0;
        end
        
        function x_ini = compute_initial_value(obj,x_ini,cost,constraint)
            % !! REMOVE WHEN DesginVariableInitializer CONSIDERS
            % Projected_Slerp INITIAL GUESS !!
        end
        
    end
    
    methods (Access = protected)
        
        function print(obj,x,iter,cost,constraint)
            if (obj.printing)
                obj.printer.print(x,iter,cost,constraint);
            end
        end
        
        function writeToFile(obj,nstep,cost,constraint)
            if ~(obj.printing)
                return
            end
            if obj.niter==1
                msh_file = fullfile('Output',obj.fileName,strcat(obj.fileName,'.txt'));
                fid_mesh = fopen(msh_file,'wt');
            else
                msh_file = fullfile('Output',obj.fileName,strcat(obj.fileName,'.txt'));
                fid_mesh = fopen(msh_file,'at');
            end
            fprintf(fid_mesh,'-----------------------------------------------------------------------------------------------\n');
            fprintf(fid_mesh,'\n');
            fprintf(fid_mesh,'Iteration: %i \n',obj.niter);
            fprintf(fid_mesh,'Nstep: %i \n',nstep);
            fprintf(fid_mesh,'Cost  %f \n',cost.value);
            for i=1:length(cost.ShapeFuncs)
                fprintf(fid_mesh,strcat('-Cost ',num2str(i),': %f \n'),cost.ShapeFuncs{i}.value);
            end
            fprintf(fid_mesh,'Constraint: %f \n',constraint.value);
            for i=1:length(constraint.ShapeFuncs)
                fprintf(fid_mesh,strcat('-Constraint ',num2str(i),': %f \n'),constraint.ShapeFuncs{i}.value);
            end
            
            switch obj.optimizer
                case {'SLERP','PROJECTED GRADIENT','HAMILTON-JACOBI'}
                    fprintf(fid_mesh,'Optimality tolerance: %f \n',obj.optimizer_unconstr.opt_cond);
                    fprintf(fid_mesh,'Kappa: %f \n',obj.optimizer_unconstr.line_search.kappa);
                case 'MMA'
                    fprintf(fid_mesh,'Optimality tolerance: %f \n',obj.kktnorm);
                case 'IPOPT'
                    fprintf(fid_mesh,'Optimality tolerance: %f \n',obj.data.inf_du);
            end
            fprintf(fid_mesh,'\n');
            fclose(fid_mesh);
        end
        
        function createOptimizerPrinter(obj,cost,constraint)
            m  = obj.mesh;            
            fN = obj.fileName;
            op = obj.optimizer;
            pM = obj.printMode;
            obj.printer = OptimizerPrinter.create(m,op,fN,pM,cost,constraint);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,settings,mesh)
            obj.fileName    = settings.case_file;
            obj.maxiter     = settings.maxiter;
            obj.printing    = settings.printing;
            obj.optimizer   = settings.optimizer;
            obj.printMode = settings.printMode;
            obj.mesh     = mesh;
        end
               
        function printFinal(obj,x,cost,constraint)
            if obj.monitoring.plotting_ON
                if obj.printing
                    obj.print(x,obj.niter,cost,constraint);
                else
                    obj.printing = 1;
                    obj.print(x,obj.niter,cost,constraint);
                    obj.printing = 0;
                end
            end            
        end
        
    end
    
end
