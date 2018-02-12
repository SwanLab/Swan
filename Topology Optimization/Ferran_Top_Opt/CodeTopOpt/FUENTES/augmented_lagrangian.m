function [x,fval,lambda,penalty,iter_level_set] = augmented_lagrangian(fobj,x0,constr_fun,options,varargin)

% AUGMENTED LAGRANGIAN OPTIMIZER
% Uses projected gradient and level set design variable update algorithms
% 
% Input:
%   fobj: Objective function handle -> [fun,gradient_fun] = fobj(x)
%   x0: Initial guess
%   constr_fun: Constraint/s function handle (inequality) -> [c,[],gradient_c,[]] = constr_fun(x)
%   options
%       - lambda: Lagrange multipliers of the constraints (array of same length as constr_fun)
%       - penalty: Penalty multipliers of the constraints (array of same length as constr_fun)
%       - update_penalty: Function handle for updating the penalty
%       - OutputFcn: Function handle to execute after each iteration (optional)
%       - optimality_tol: Stopping optimality criteria (theta for level set and incre_gamma for projected gradient)
%       - constr_tol: Tolerance for the constraints, can be a scalar or an array of same length as constr_fun
%       - algorithm: 'Projected_gradient' or 'level_set'
%       - kappa_min: Minimum kappa to use
%       - scalar_product: Function handle for computing scalar products
%       - kappa0: Initial kappa value (optional)
%       - lb: Lower bound for the design variable (only projected gradient)
%       - ub: Upper bound for the design variable (only projected gradient)
%       - constraint_case: equality (c = 0) or inequality (c < 0)
%   varargin: Additional data to be passed to OutputFcn
% 
% Example
% [x,fval,lambda,penalty] = augmented_lagrangian(fobj,x0,constr_fun,options,varargin)

%% Initialize problem and options
design_variable = x0;

options.call_theta = @(phi,g) call_theta_level_set(phi,g,options.scalar_product,options.algorithm);
options.update_lambda = @(lambda,penalty,constraint) lambda + penalty.*constraint';
lambda = options.lambda(:)'; % horizontal vector
penalty = options.penalty(:)'; % horizontal vector

nconstr = length(lambda);
nconstr_tol = length(options.constr_tol);
if nconstr_tol < nconstr
   options.constr_tol(nconstr_tol+1:nconstr) = 0; % default
end

%% Start algorithm
[compliance,gradient_compliance] = fobj(design_variable);
[constraints,~,gradient_constraints,~] = constr_fun(design_variable);
[constraints,gradient_constraints] = check_constraint_case(constraints,gradient_constraints,options.constraint_case,lambda,penalty);

cost = assemble_cost(compliance,constraints,lambda,penalty);
gradient = assemble_gradient(gradient_compliance,constraints,gradient_constraints,lambda,penalty);
theta = options.call_theta(design_variable,gradient); 

incre_x = 1;
iter_level_set = 0;
kappa = 0;
design_variable_old = design_variable;

if ~isempty(options.OutputFcn) % call output functions (update post_info, plots and save to file)
    optdata.cost = cost;
    optdata.theta = theta;
    optdata.lambda = lambda;
    optdata.penalty = penalty;
    optdata.kappa = kappa;
    optdata.gradient = gradient;
    optdata.iter_ls = iter_level_set;
    optdata.incre_gamma = incre_x;
    for i = 1:length(options.OutputFcn)
        feval(options.OutputFcn{i},design_variable,optdata,varargin{:});
    end
end

% stopping_criteria_opt = compute_stopping_criteria_new(theta,incre_x,constraints,options.optimality_tol,options.constr_tol,options.algorithm,penalty);
if ~isfield(options,'kappa0')
    kappa = [];
else
    kappa = options.kappa0;
end
stopping_criteria_opt = true;
while stopping_criteria_opt
    iter_level_set = iter_level_set + 1;
    
    % Maximization thanks to lambda
    penalty = options.update_penalty(penalty);
    lambda = options.update_lambda(lambda,penalty,constraints);

    cost = assemble_cost(compliance,constraints,lambda,penalty);
    gradient = assemble_gradient(gradient_compliance,constraints,gradient_constraints,lambda,penalty);
    theta = options.call_theta(design_variable,gradient); 
    
    % Minimization thanks to design_variable (Line search)
    stopping_Criteria_ls = 1; kfrac = 2;
    kappa = inicialize_kappa(design_variable,gradient,kappa*kfrac,options.scalar_product,options.algorithm);
    volume = options.volume(design_variable);
    while stopping_Criteria_ls
        design_variable_ls = design_variable_update(design_variable,gradient,theta,kappa,options);
        
        % Update augmented Lagrangian function        
        [compliance_ls,gradient_compliance_ls] = fobj(design_variable_ls);
        [constraints_ls,~,gradient_constraints_ls,~] = constr_fun(design_variable_ls);
        [constraints_ls,gradient_constraints_ls] = check_constraint_case(constraints_ls,gradient_constraints_ls,options.constraint_case,lambda,penalty);

        cost_ls = assemble_cost(compliance_ls,constraints_ls,lambda,penalty);
        gradient_ls = assemble_gradient(gradient_compliance_ls,constraints_ls,gradient_constraints_ls,lambda,penalty);
        theta_ls = options.call_theta(design_variable_ls,gradient_ls);
        volume_ls = options.volume(design_variable_ls);
        incr_vol_ls = abs(volume_ls - volume);
        %volume = volume_ls;

        kappa = kappa/kfrac;
        incr_cost = (cost_ls - cost)/abs(cost);
%         stopping_Criteria_ls = ~(incr_cost < 0 || kappa <= options.kappa_min);
        stopping_Criteria_ls = ~((incr_cost < 0 && incr_vol_ls < options.max_constr_change) || kappa <= options.kappa_min);
    end

    
    if kappa > options.kappa_min
        compliance = compliance_ls;
        gradient_compliance = gradient_compliance_ls;
        design_variable = design_variable_ls; 
        cost = cost_ls; 
        gradient = gradient_ls;
        constraints = constraints_ls;
        gradient_constraints = gradient_constraints_ls;
        theta = theta_ls;
    end

    switch options.algorithm
    case 'level_set'
        incre_x = 1;
    case 'Projected_gradient'
        incre_x = sqrt(options.scalar_product(design_variable - design_variable_old,design_variable - design_variable_old))/sqrt(options.scalar_product(design_variable_old,design_variable_old));
    end

    if ~isempty(options.OutputFcn) % call output functions (update post_info, plots and save to file)
        optdata.cost = cost;
        optdata.theta = theta;
        optdata.lambda = lambda;
        optdata.penalty = penalty;
        optdata.kappa = kappa*kfrac;
        optdata.compliance = compliance;
        optdata.gradient = gradient;
        optdata.iter_ls = iter_level_set;
        optdata.incre_gamma = incre_x;
        for i = 1:length(options.OutputFcn)
            feval(options.OutputFcn{i},design_variable,optdata,varargin{:});
        end
    end
    
    stopping_criteria_opt = compute_stopping_criteria_new(theta,incre_x,constraints,options.optimality_tol,options.constr_tol,options.algorithm,penalty);
    design_variable_old = design_variable;

end

x = design_variable;
fval = compliance;

end


function  cost = assemble_cost(compliance,constraints,lambda,penalty)

cost = compliance + lambda*constraints + 0.5*penalty*(constraints.*constraints);

end


function gradient = assemble_gradient(gradient_compliance,constraints,gradient_constraints,lambda,penalty)

% Note: the gradient of each constraint must be in column form, gradient
% constraints is a matrix of nvariables x nconstr

gradient = gradient_constraints*lambda' + gradient_constraints*(penalty'.*constraints) + gradient_compliance;

end

function [constraints,gradient_constraints] = check_constraint_case(constraints0,gradient_constraints0,constraint_case,lambda,penalty)

switch constraint_case
    case 'equality'
        constraints = constraints0;
        gradient_constraints = gradient_constraints0;
        
    case 'inequality'
        contr_inactive_value = -lambda(:)./penalty(:);
        inactive_constr = contr_inactive_value > constraints0;
        constraints = constraints0;
        constraints(inactive_constr) = contr_inactive_value(inactive_constr);
        gradient_constraints = gradient_constraints0;
        gradient_constraints(:,inactive_constr) = 0;
        
    otherwise
        error('Constraint case not valid.');
end



end