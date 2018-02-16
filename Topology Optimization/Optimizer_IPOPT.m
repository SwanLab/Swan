classdef Optimizer_IPOPT < Optimizer
    properties
        m
        info
        max_iter
        constraint_tolerance
        optimality_tolerance
    end 
    methods
        function obj=Optimizer_IPOPT(settings)
            obj@Optimizer(settings);
            obj.m=settings.nconstr;
            obj.max_iter=5e3;
        end
        function optimality_tolerance=get.optimality_tolerance(obj)
            optimality_tolerance=obj.target_parameters.optimality_tol;
        end
        function constraint_tolerance=get.constraint_tolerance(obj)
            constraint_tolerance=obj.target_parameters.constr_tol*1e-1;
        end
        
        function x = solveProblem(obj,x_ini,cost,constraint,interpolation,filter) 
            cost.computef(x_ini,obj.physicalProblem,interpolation,filter)
            funcs.objective = @(x) obj.objective(x,cost,interpolation,filter);
            funcs.gradient = @(x) obj.gradient(x,cost,obj.physicalProblem,interpolation,filter);
            funcs.constraints = @(x) obj.constraint(x,constraint,obj.physicalProblem,interpolation,filter);
            funcs.jacobian = @(x) sparse(obj.constraint_gradient(x,constraint,obj.physicalProblem,interpolation,filter)');
            n = length(x_ini);
            funcs.jacobianstructure = @() sparse(ones(obj.m,n));
            plotx=@(x) obj.plotX(x);
            funcs.iterfunc = @(iter,fval,data) obj.outputfun_ipopt(iter,fval,data,plotx);
            
            options.ipopt.print_level           = 0;
            options.ipopt.hessian_approximation = 'limited-memory';
            options.ipopt.limited_memory_update_type = 'bfgs';
            options.ub = ones(length(x_ini),1);
            options.lb = zeros(length(x_ini),1);
            if isfield(options,'constraint_case')
                if strcmp(options.constraint_case,'equality')
                    options.cl = zeros(obj.m,1);
                else
                    options.cl = -Inf*ones(obj.m,1); % lower bound constraint
                end
            else
                options.cl = -Inf*ones(obj.m,1);
            end
            options.cu = zeros(obj.m,1); % upper bound constraint value
            options.ipopt.max_iter = obj.max_iter;
            options.ipopt.constr_viol_tol = obj.constraint_tolerance;
            %         options.ipopt.dual_inf_tol = optimality_tolerance;
            options.ipopt.compl_inf_tol = obj.constraint_tolerance;
            options.ipopt.tol = obj.optimality_tolerance;
            
            [x, obj.info] = ipopt(x_ini,funcs,options);
        end
        function f=objective(obj,x,cost,interpolation,filter)
            cost.computef(x,obj.physicalProblem,interpolation,filter)
            f=cost.value;
        end
    end
    methods (Static)
        function f=constraint(x,constraint,physProblem,interpolation,filter)
            constraint.computef(x,physProblem,interpolation,filter)
            f=constraint.value;
        end
        function g=gradient(x,cost,physProblem,interpolation,filter)
            cost.computef(x,physProblem,interpolation,filter)
            g=cost.gradient;
        end
        function g=constraint_gradient(x,constraint,physProblem,interpolation,filter)
            constraint.computef(x,physProblem,interpolation,filter)
            g=constraint.gradient;
        end
        function stop=outputfun_ipopt(iter,fval,data,plotx)
            disp(strcat('Iter:',num2str(iter)));
            stop=true;
            plotx(data.x);
        end
    end
end