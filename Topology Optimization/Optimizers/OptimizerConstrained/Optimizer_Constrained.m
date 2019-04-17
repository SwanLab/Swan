classdef Optimizer_Constrained < Optimizer
    
    properties (Access = public)
        fhtri
        niter = 0
        optimizer
        maxiter
    end
    
    properties (Access = protected)
        monitor
        cost
        constraint
        constraintCase
    end
    
    properties (Access = private)
        designVar
        postProcess
        printing
        printMode
    end
    
    methods (Access = public)
        
        function obj = Optimizer_Constrained(cParams)
            set = cParams.settings.settings;
            designVar = cParams.designVariable;
            
            obj.constraintCase   = set.constraint_case;
            obj.hasConverged      = false;
            
            % !! REMOVE !!
            obj.cost = cParams.settings.cost;
            obj.constraint = cParams.settings.constraint;
            
            obj.init(set,designVar);
            
            mS.showOptParams = cParams.monitoring;
            mS.plotting  = set.plotting;
            mS.settings  = set;
            mS.designVar = cParams.designVariable;
            
            obj.monitor = MonitoringDocker(mS);
        end
        
        function designVar = solveProblem(obj,designVar,cost,constraint,istep,nstep)
            obj.createPostProcess(cost,constraint);
            x_ini = designVar.value;
            cost.computeCostAndGradient(x_ini);
            constraint.computeCostAndGradient(x_ini);
            obj.print(x_ini,obj.niter,cost,constraint);
            
            obj.monitor.refresh(x_ini,obj.niter,cost,constraint,obj.stop_vars,obj.hasFinished(istep,nstep),istep,nstep);
            
            while ~obj.hasFinished(istep,nstep)
                obj.niter = obj.niter+1;
                x = obj.update(x_ini);
                obj.monitor.refresh(x,obj.niter,cost,constraint,obj.stop_vars,obj.hasFinished(istep,nstep),istep,nstep);
                obj.print(x,obj.niter,cost,constraint);
                obj.writeToFile(istep,cost,constraint)
                x_ini = x;
            end
            obj.printFinal(x,cost,constraint);
            designVar.update(x);
            
            obj.hasConverged = 0;
        end
        
    end
    
    methods (Access = protected)
        
        function print(obj,x,iter,cost,constraint)
            if (obj.printing)
                d.x = x;
                d.cost = cost;
                d.constraint = constraint;
                obj.postProcess.print(iter,d);
            end
        end
        
        function writeToFile(obj,nstep,cost,constraint)
            if ~(obj.printing)
                return
            end
            fileName = obj.designVar.meshGiD.problemID;
            if obj.niter == 1
                msh_file = fullfile('Output',fileName,strcat(fileName,'.txt'));
                fid_mesh = fopen(msh_file,'wt');
            else
                msh_file = fullfile('Output',fileName,strcat(fileName,'.txt'));
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
        
        function createPostProcess(obj,cost,constraint)
            fileName = obj.designVar.meshGiD.problemID;
            d = obj.createPostProcessDataBase(fileName);
            d.printMode = obj.printMode;
            d.optimizer = obj.optimizer;
            d.cost = cost;
            d.constraint = constraint;
            obj.postProcess = Postprocess('TopOptProblem',d);
        end
        
    end
    
    methods (Access = protected)
        
        function itHas = hasFinished(obj,istep,nstep)
            itHas = obj.hasConverged || obj.niter >= obj.maxiter*(istep/nstep);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,settings,designVar)
            obj.designVar = designVar;
            obj.optimizer = settings.optimizer;
            obj.maxiter   = settings.maxiter;
            obj.printing  = settings.printing;
            obj.printMode = settings.printMode;
        end
        
        function d = createPostProcessDataBase(obj,fileName)
            d.mesh    = obj.designVar.meshGiD;
            d.outName = fileName;
            ps = PostProcessDataBaseCreator(d);
            d = ps.getValue();
        end
        
        function printFinal(obj,x,cost,constraint)
            if obj.monitor.shallDisplayDesignVar
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
