function [ fobj,constr_fun,problembsc,initial_values,options ] = create_optimization_problem( problem,problem_type,problembsc,target_parameters,initial_values,options )

global post_info
Vtarget = target_parameters.Vtarget;
epsilon = target_parameters.epsilon;
epsilon_isotropy = target_parameters.epsilon_isotropy;
Cstar = target_parameters.Cstar;

% Penalization parameter (perimeter)
lambda = 0.0;
    
% First time
if isempty(initial_values)
    % Input data for MICRO
    problembsc.costfunct = 'alpha'; % 'alpha' 'fraction' 'C-C*'
    problembsc.enforceCh_type = 'CCstar'; % 'isotropy' 'CCstar'
    problembsc.alpha = [1,1,0];
    problembsc.beta = [1,1,0];
    problembsc.beta = problembsc.beta(:)/norm(problembsc.beta(:));
    problembsc.alpha = problembsc.alpha(:)/norm(problembsc.alpha(:));
    
    % Cstar
    problembsc.Ch_star_final = compute_Ch_star(problem.matprop);
    problembsc.Ch_star = problembsc.Ch_star_final;
    problembsc.selectiveC_Cstar = [1,1,1;
                                   1,1,1;
                                   1,1,1]; % 1 to select the component
    
    % Adimensional values
    compliancefun = @(x) problem.compliance(x,problembsc,1);
    [initial_values.compliance0,~] = compliancefun(problem.x0);
    structural_values = post_info.structural_values;
    %[h_C_0,~,~] = problem.equilibrium(problem.gamma_hC0,problembsc,1);
    %[~,~,structural_values] = problem.equilibrium(problem.x0_elemental,problembsc,1);
    Cstar = structural_values.matCh;
    post_info.Ch_star = problembsc.Ch_star;
    enforceChdiff = @(x,epsilon_isotropy,component) problem.enforceChdiff_square(x,epsilon_isotropy,compliancefun,problembsc.enforceCh_type,component,ones(6,1));
    Lsquare_fun = @(x) Lsquare_Ch_difference (x,enforceChdiff,0,1);
    if strcmp(problem.design_type,'MICRO')
        initial_values.Lsquare0 = Lsquare_fun(problem.x0);
        initial_values.componentCh = ones(6,1);
    else
        initial_values.Lsquare0 = 1;
        initial_values.componentCh = ones(6,1);
    end
    
    % Initialize options
    options = struct;
end
problembsc.Ch_star = Cstar;
post_info.Ch_star = problembsc.Ch_star;

% Define problem functions
switch problem_type
    case 'min_compliance_st_vol_per'
        fobj = @(x) problem.compliance(x,problembsc,initial_values.compliance0);
        constr_cell{1} = @(x) problem.volume(x,Vtarget);
        constr_cell{2} = @(x) problem.perimeter(x,epsilon);
        constr_fun = @(x) constr_opt_fun(x,constr_cell(:));
        
        % Post-process
        post_info.plot_names = {'Volume','Perimeter'}; % constraint names
        
    case 'min_compliance+lam*per_st_vol'      
        perimeterfun = @(x) problem.perimeter(x,epsilon);
        fobj = @(x) compliance_perimeter( x,@(x) problem.compliance(x,problembsc,initial_values.compliance0),perimeterfun,lambda );
        constr_cell{1} = @(x) problem.volume(x,Vtarget);
        constr_fun = @(x) constr_opt_fun(x,constr_cell(:));
        
        % Post-process
        post_info.plot_names = {'Volume','Perimeter'}; % constraint names
        post_info.lambda_fobj = lambda;
        
    case 'min_weight_st_C-C*'       
        fobj = @(x) problem.volume(x,Vtarget);
        compliancefun = @(x) problem.compliance(x,problembsc,initial_values.compliance0);
        enforceChdiff = @(x,epsilon_isotropy,component) problem.enforceChdiff_square(x,0,compliancefun,problembsc.enforceCh_type,component,ones(6,1));
        Lsquare_fun = @(x) Lsquare_Ch_difference (x,enforceChdiff,epsilon_isotropy,initial_values.Lsquare0);
        constr_cell{1} = Lsquare_fun;
        constr_fun = @(x) constr_opt_fun(x,constr_cell(:));
        
        % Initialize lambda for C-C*
        if ~isfield(options,'lambda')
            m = numel(constr_cell);
            options.lambda = zeros(1,m);
            options.lambda(1) = 2000;
        end
        
        % Post-process
        post_info.plot_names = {'||C - Cstar||2','Perimeter'}; % constraint names
        
    case 'min_compliance+lam*per_st_vol_enforceCh_inf'      
        perimeterfun = @(x) problem.perimeter(x,epsilon);
        compliancefun = @(x) problem.compliance(x,problembsc,initial_values.compliance0);
        fobj = @(x) compliance_perimeter( x,compliancefun,perimeterfun,lambda );
        constr_cell{1} = @(x) problem.volume(x,Vtarget);
        constr_cell{2} = @(x) problem.enforceChdiff_square(x,epsilon_isotropy,compliancefun,problembsc.enforceCh_type,1,initial_values.componentCh);
        constr_cell{3} = @(x) problem.enforceChdiff_square(x,epsilon_isotropy,compliancefun,problembsc.enforceCh_type,2,initial_values.componentCh);
        constr_cell{4} = @(x) problem.enforceChdiff_square(x,epsilon_isotropy,compliancefun,problembsc.enforceCh_type,3,initial_values.componentCh);
        constr_cell{5} = @(x) problem.enforceChdiff_square(x,0,compliancefun,problembsc.enforceCh_type,4,initial_values.componentCh);
        constr_cell{6} = @(x) problem.enforceChdiff_square(x,0,compliancefun,problembsc.enforceCh_type,5,initial_values.componentCh);
        constr_cell{7} = @(x) problem.enforceChdiff_square(x,epsilon_isotropy,compliancefun,problembsc.enforceCh_type,6,initial_values.componentCh);
        constr_fun = @(x) constr_opt_fun(x,constr_cell(:));
        
        % Post-process
        post_info.plot_names = {'Volume and Perimeter',problembsc.enforceCh_type}; % constraint names
        post_info.constr_names1 = {'Volume','Perimeter'};
        post_info.constr_names2 = cellstr(('1':'6')');
        post_info.lambda_fobj = lambda;
        
    case 'min_compliance+lam*per_st_enforceCh_inf_equality'      
        perimeterfun = @(x) problem.perimeter(x,epsilon);
        compliancefun = @(x) problem.compliance(x,problembsc,initial_values.compliance0);
        fobj = @(x) compliance_perimeter( x,compliancefun,perimeterfun,lambda );
        constr_cell{1} = @(x) problem.enforceChdiff(x,compliancefun,problembsc.enforceCh_type,1,initial_values.componentCh);
        constr_cell{2} = @(x) problem.enforceChdiff(x,compliancefun,problembsc.enforceCh_type,2,initial_values.componentCh);
        constr_cell{3} = @(x) problem.enforceChdiff(x,compliancefun,problembsc.enforceCh_type,3,initial_values.componentCh);
        constr_cell{4} = @(x) problem.enforceChdiff(x,compliancefun,problembsc.enforceCh_type,4,initial_values.componentCh);
        constr_cell{5} = @(x) problem.enforceChdiff(x,compliancefun,problembsc.enforceCh_type,5,initial_values.componentCh);
        constr_cell{6} = @(x) problem.enforceChdiff(x,compliancefun,problembsc.enforceCh_type,6,initial_values.componentCh);
        constr_fun = @(x) constr_opt_fun(x,constr_cell(:));
        
      % Initialize lambda
        if ~isfield(options,'lambda')
            m = numel(constr_cell);
            options.lambda = 0*ones(1,m);
        end
        if ~isfield(options,'constraint_case')
            options.constraint_case = 'equality'; % 'equality' 'inequality'
        end
                
        % Post-process
        post_info.plot_names = {'Volume and Perimeter',problembsc.enforceCh_type}; % constraint names
        post_info.constr_names1 = {'Volume','Perimeter'};
        post_info.constr_names2 = cellstr(('1':'6')');
        post_info.lambda_fobj = lambda;
        
    case 'min_compliance+lam*per_st_enforceCh_equality'      
        perimeterfun = @(x) problem.perimeter(x,epsilon);
        compliancefun = @(x) problem.compliance(x,problembsc,initial_values.compliance0);
        fobj = @(x) compliance_perimeter( x,compliancefun,perimeterfun,lambda );
        enforceChdiff = @(x,component) problem.enforceChdiff(x,compliancefun,problembsc.enforceCh_type,component,ones(6,1));
        constr_cell{1} = @(x) Lequality_Ch_difference (x,enforceChdiff,1);
        constr_fun = @(x) constr_opt_fun(x,constr_cell(:));
        
      % Initialize lambda
        if ~isfield(options,'lambda')
            m = numel(constr_cell);
            options.lambda = 0*ones(1,m);
        end
        if ~isfield(options,'constraint_case')
            options.constraint_case = 'equality'; % 'equality' 'inequality'
        end
                
        % Post-process
        post_info.plot_names = {'Volume and Perimeter',problembsc.enforceCh_type}; % constraint names
        post_info.constr_names1 = {'Volume','Perimeter'};
        post_info.constr_names2 = cellstr(('1':'6')');
        post_info.lambda_fobj = lambda;
        
    case 'min_compliance+lam*per_st_enforceCh_L2'      
        perimeterfun = @(x) problem.perimeter(x,epsilon);
        compliancefun = @(x) problem.compliance(x,problembsc,initial_values.compliance0);
        fobj = @(x) compliance_perimeter( x,compliancefun,perimeterfun,lambda );
        enforceChdiff = @(x,epsilon_isotropy,component) problem.enforceChdiff_square(x,0,compliancefun,problembsc.enforceCh_type,component,ones(6,1));
        constr_cell{1} = @(x) Lsquare_Ch_difference (x,enforceChdiff,epsilon_isotropy,1);
        constr_fun = @(x) constr_opt_fun(x,constr_cell(:));
                
        % Post-process
        post_info.plot_names = {problembsc.enforceCh_type,' '}; % constraint names
        post_info.lambda_fobj = lambda;
        
    case 'min_weight_st_enforceCh_inf'
        volumefun = @(x) problem.volume(x,Vtarget);
        perimeterfun = @(x) problem.perimeter(x,epsilon);
        fobj = @(x) compliance_perimeter( x,volumefun,perimeterfun,lambda );
        compliancefun = @(x) problem.compliance(x,problembsc,initial_values.compliance0);
        constr_cell{1} = @(x) problem.enforceChdiff_square(x,epsilon_isotropy,compliancefun,problembsc.enforceCh_type,1,initial_values.componentCh);
        constr_cell{2} = @(x) problem.enforceChdiff_square(x,epsilon_isotropy,compliancefun,problembsc.enforceCh_type,2,initial_values.componentCh);
        constr_cell{3} = @(x) problem.enforceChdiff_square(x,epsilon_isotropy,compliancefun,problembsc.enforceCh_type,3,initial_values.componentCh);
        constr_cell{4} = @(x) problem.enforceChdiff_square(x,0,compliancefun,problembsc.enforceCh_type,4,initial_values.componentCh);
        constr_cell{5} = @(x) problem.enforceChdiff_square(x,0,compliancefun,problembsc.enforceCh_type,5,initial_values.componentCh);
        constr_cell{6} = @(x) problem.enforceChdiff_square(x,epsilon_isotropy,compliancefun,problembsc.enforceCh_type,6,initial_values.componentCh);
        constr_fun = @(x) constr_opt_fun(x,constr_cell(:));
        
        % Initialize lambda for C-C*
        if ~isfield(options,'lambda')
            m = numel(constr_cell);
            options.lambda = 0*ones(1,m);
        end
        
        % Post-process
        post_info.plot_names = {'Perimeter',problembsc.enforceCh_type}; % constraint names
        post_info.constr_names1 = {'Perimeter'};
        post_info.constr_names2 = {'C11','C22','C33','C23','C13','C12'};
        post_info.lambda_fobj = lambda;
        
    case 'min_enforceCh_L2+lam*per_st_vol'      
        perimeterfun = @(x) problem.perimeter(x,epsilon);
        compliancefun = @(x) problem.compliance(x,problembsc,initial_values.compliance0);
        enforceChdiff = @(x,epsilon_isotropy,component) problem.enforceChdiff_square(x,0,compliancefun,problembsc.enforceCh_type,component,ones(6,1));
        Lsquare_fun = @(x) Lsquare_Ch_difference (x,enforceChdiff,epsilon_isotropy,initial_values.Lsquare0);
        fobj = @(x) compliance_perimeter( x,Lsquare_fun,perimeterfun,lambda );
        constr_cell{1} = @(x) problem.volume(x,Vtarget);
        constr_fun = @(x) constr_opt_fun(x,constr_cell(:));
        
        % Post-process
        post_info.plot_names = {'Volume',problembsc.enforceCh_type}; % constraint names
        post_info.lambda_fobj = lambda;
          
        
    case 'min_weight_st_enforceCh_inf_equality'
        volumefun = @(x) problem.volume(x,Vtarget);
        perimeterfun = @(x) problem.perimeter(x,epsilon);
        fobj = @(x) compliance_perimeter( x,volumefun,perimeterfun,lambda );
        compliancefun = @(x) problem.compliance(x,problembsc,initial_values.compliance0);
        constr_cell{1} = @(x) problem.enforceChdiff(x,compliancefun,problembsc.enforceCh_type,1,initial_values.componentCh);
        constr_cell{2} = @(x) problem.enforceChdiff(x,compliancefun,problembsc.enforceCh_type,2,initial_values.componentCh);
        constr_cell{3} = @(x) problem.enforceChdiff(x,compliancefun,problembsc.enforceCh_type,3,initial_values.componentCh);
        constr_cell{4} = @(x) problem.enforceChdiff(x,compliancefun,problembsc.enforceCh_type,4,initial_values.componentCh);
        constr_cell{5} = @(x) problem.enforceChdiff(x,compliancefun,problembsc.enforceCh_type,5,initial_values.componentCh);
        constr_cell{6} = @(x) problem.enforceChdiff(x,compliancefun,problembsc.enforceCh_type,6,initial_values.componentCh);
        constr_fun = @(x) constr_opt_fun(x,constr_cell(:));
        
        % Initialize lambda for C-C*
        if ~isfield(options,'lambda')
            m = numel(constr_cell);
            options.lambda = 0*ones(1,m);
        end
        if ~isfield(options,'constraint_case')
            options.constraint_case = 'equality'; % 'equality' 'inequality'
        end
        
        % Post-process
        post_info.plot_names = {'Perimeter',problembsc.enforceCh_type}; % constraint names
        post_info.constr_names1 = {'Perimeter'};
        post_info.constr_names2 = {'C11','C22','C33','C23','C13','C12'};
        post_info.lambda_fobj = lambda;
end
problembsc.nconstr = numel(constr_cell);

end

