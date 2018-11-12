function [design_variable,fval,options] = optimizer_selection (fobj,constr_fun,x0,problembsc,problem,h_C_0,options,optimality_tolerance,constraint_tolerance)
global post_info

% General options
m = problembsc.nconstr; % number of constraints
max_iter = 5000;

switch problem.algorithm
    case {'level_set','Projected_gradient'} % use augmented Lagrangian
        % Options
        if ~isfield(options,'lambda')
            options.lambda = 0*ones(1,m);
        end
        if ~isfield(options,'penalty')
            options.penalty = 1*ones(1,m);
        end
        options.update_penalty = @(penalty) penalty;%min(10,(1.15)^(1)*penalty);
        if strcmp(problem.algorithm,'level_set')
            options.optimality_tol = optimality_tolerance*pi/180;
        else
            options.optimality_tol = optimality_tolerance;
        end
        options.constr_tol = constraint_tolerance*ones(m,1);
        options.algorithm = problem.algorithm;
        options.kappa_min = 1e-15;
        options.max_constr_change = +Inf;
        options.scalar_product = problem.scalar_product;
        options.lb = problem.lb;
        options.ub = problem.ub;
        if ~isfield(options,'constraint_case')
            options.constraint_case = 'inequality'; % 'equality' 'inequality'
        end
        options.volume = problem.volume_current;
        options.OutputFcn = {problem.outputopt};
        
        % Algorithm
        [design_variable,fval,options.lambda,options.penalty,global_iterations] = augmented_lagrangian(fobj,x0,constr_fun,options,false);
        post_info.cumulative_iter = post_info.cumulative_iter + global_iterations;
        
    case {'MMA','GCMMA'}
        % Options
        ub = problem.ub;
        lb = problem.lb;
        options.outeriter = 0;
        options.maxoutit = max_iter;
        options.kkttol = optimality_tolerance;
        if ~isfield(options,'constraint_case')
            options.constraint_case = 'inequality'; % 'equality' 'inequality'
        end
        if strcmp(options.constraint_case,'equality')
            m = 2*m;
        end
        c_mma = 1e3*ones(m,1);
        d_mma = 0*ones(m,1);
        a0 = 1;
        a = 0*ones(m,1);
        options.OutputFcn = {problem.outputopt};
        
        % Algorithm
        fun = @(x) funmma(x,fobj,constr_fun,options.constraint_case);
        switch problem.algorithm
            case 'MMA'
                [design_variable,fval,constr_val,kktnorm,~,outit] = mma(fun,x0,lb,ub,x0,c_mma,d_mma,a0,a,options,false);
            case 'GCMMA'
                [design_variable,fval,constr_val,kktnorm,~,outit] = gcmma(fun,x0,lb,ub,x0,c_mma,d_mma,a0,a,options,false);
        end
        post_info.cumulative_iter = post_info.cumulative_iter + outit;
        
    case 'fmincon'
        % Options
        lb = problem.lb;
        ub = problem.ub;
        options = optimoptions('fmincon', ...
                               options, ...
                               'Algorithm','interior-point',...
                               'ConstraintTolerance',constraint_tolerance*ones(m,1), ...
                               'DiffMaxChange',Inf, ... % maximum change in variables
                               'Display','none', ...
                               'FunctionTolerance',1e-3, ...
                               'FunValCheck','off', ... % check if fobj values are complex, Inf or NaN
                               'MaxFunctionEvaluations',4*max_iter, ...
                               'MaxIterations',max_iter, ...
                               'OptimalityTolerance',optimality_tolerance, ...
                               'OutputFcn',{problem.outputopt}, ...
                               'SpecifyConstraintGradient',true, ...
                               'SpecifyObjectiveGradient',true, ...
                               'StepTolerance',1e-3); % incre gamma
                           
        switch options.Algorithm
            case 'interior-point'
                options = optimoptions(options, ... % update current options for specific algorithm
                                       'HessianApproximation','bfgs', ...
                                       'HonorBounds',true, ... % bound constraints are satisfied at every iteration
                                       'InitBarrierParam',0.1, ... % larger values might be helpful
                                       'InitTrustRegionRadius',sqrt(length(x0)), ... % sometimes smaller values might help
                                       'ScaleProblem','none', ...
                                       'SubproblemAlgorithm','factorization');
                                   
            otherwise
                error('fmincon''s algorithm %s not detected!',options.Algorithm);
        end
        
        % Algorithm
        [design_variable,fval,exitflag,output,lambda,grad,hessian] = fmincon(fobj,x0,[],[],[],[],lb,ub,constr_fun,options);
        post_info.cumulative_iter = post_info.cumulative_iter + output.iterations;
    
    case 'IPOPT'
        post_info.xold_fobj = zeros(size(x0));
        post_info.xold_constr = zeros(size(x0));
        funcs.objective = @(x) fobj_ipopt(x,fobj);
        funcs.gradient = @(x) fobj_gradient_ipopt(x,fobj);
        funcs.constraints = @(x) constraints_ipopt(x,constr_fun);
        funcs.jacobian = @(x) constraints_gradient_ipopt(x,constr_fun);
        n = length(x0);
        funcs.jacobianstructure = @() sparse(ones(m,n));
        funcs.iterfunc = @(iter,fval,data) outputfun_ipopt(iter,fval,data,problem.outputopt);
        
        % Options
        options.ipopt.print_level           = 0;
        options.ipopt.hessian_approximation = 'limited-memory';
        options.ipopt.limited_memory_update_type = 'bfgs';
        options.ub = problem.ub;
        options.lb = problem.lb;
        if isfield(options,'constraint_case')
            if strcmp(options.constraint_case,'equality')
                options.cl = zeros(m,1);
            else
                options.cl = -Inf*ones(m,1); % lower bound constraint
            end
        else
            options.cl = -Inf*ones(m,1);
        end
        options.cu = zeros(m,1); % upper bound constraint value
        options.ipopt.max_iter = max_iter;
        options.ipopt.constr_viol_tol = constraint_tolerance;
%         options.ipopt.dual_inf_tol = optimality_tolerance;
        options.ipopt.compl_inf_tol = constraint_tolerance;
        options.ipopt.tol = optimality_tolerance;
        
        [design_variable, info] = ipopt(x0,funcs,options);
        fval = [];
        
        % Update options for future iterations
        options.lambda = info.lambda;
        options.zl = info.zl;
        options.zu = info.zu;
        options.ipopt.warm_start_init_point = 'yes';
        
        post_info.cumulative_iter = post_info.cumulative_iter + info.iter;
        
    otherwise
        error('Algorithm not implemented.')
end

end