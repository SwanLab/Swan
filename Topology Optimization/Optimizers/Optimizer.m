classdef Optimizer < handle
    properties
        has_converged = false;
        stop_vars
        target_parameters = struct;
        nconstr
        constraint_case
        postprocess
    end
    
    properties (Access = ?Optimizer_Constrained)
        case_file
        printing
        monitoring
        mesh
    end
    
    methods
        function obj = Optimizer(settings)
            obj.nconstr = settings.nconstr;
            obj.case_file=settings.case_file;
            obj.target_parameters = settings.target_parameters;
            obj.constraint_case=settings.constraint_case;
            obj.postprocess = Postprocess_TopOpt.Create(settings.optimizer);
        end
    end
    
    methods (Abstract)
        % x = updateX(obj,x_ini,cost,constraint); %% !! IPOPT doesn't use it (black box) !!
    end
    
    methods (Access = protected)
        function cons = setConstraint_case(obj,constraint)
            cons = constraint;
            switch obj.constraint_case
                case 'EQUALITY'
                case 'INEQUALITY'
                    contr_inactive_value = -obj.objfunc.lambda(:)./obj.objfunc.penalty(:);
                    inactive_constr = contr_inactive_value' > constraint.value;
                    cons.value(inactive_constr) = contr_inactive_value(inactive_constr);
                    cons.gradient(:,inactive_constr) = 0;
                otherwise
                    error('Constraint case not valid.');
            end
        end
        function print(obj,design_variable,iter)
            if ~(obj.printing)
                return
            end
            
            %results.physicalVars = obj.physicalProblem.variables;
            results.design_variable = design_variable;
            results.iter = iter;
            results.case_file = obj.case_file;
            %results.design_variable_reg = design_variable_reg;
            obj.postprocess.print(obj.mesh,results);
        end
        function writeToFile(obj,nstep,cost,constraint)
            if ~(obj.printing)
                return
            end
            if obj.niter==1
                msh_file = fullfile('Output',obj.case_file,strcat(obj.case_file,'.txt'));
                fid_mesh = fopen(msh_file,'wt');
            else
                msh_file = fullfile('Output',obj.case_file,strcat(obj.case_file,'.txt'));
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
    end
end