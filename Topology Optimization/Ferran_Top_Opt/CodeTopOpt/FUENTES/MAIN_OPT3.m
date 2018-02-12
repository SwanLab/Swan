function  MAIN_OPT3 (file_write,file_name,Vfrac,Ptarget,method,Algorithm,kernel_case,matprop,TYPE,problem_type,incropt,init_design_variable,perimeter_case,fhplots,select_screen)

if nargin == 0
    % Optimization options   
    Vfrac=0.5;
    Ptarget=1;

    Algorithm = 'level_set';%'Projected_gradient';%'Projected_gradient';%'level_set';
    kernel_case = 'P1_kernel';%'PDE';%'P1_kernel';%'P1_kernel';'P0_kernel';
    method='SIMP_ALL'; % 'SIMP' 'SIMP_ALL'

    [file_write,file_name] = init_case_to_run;
    [fext,fext_adjoint,element,fixnodes,problembsc,coordinates,dim,Msmooth,Stiff_smooth,emass] = call_read_data(file_name,method,perimeter_case);
    [P_operator] = compute_P_operator(Msmooth,element,dim);
    regularization = 1;

    element.smoothing = 1;
else
    file_write=[file_write,'/',file_name];
    init_case_to_run(file_name);
    [fext,fext_adjoint,element,fixnodes,problembsc,coordinates,dim,Msmooth,Stiff_smooth,emass] = call_read_data(file_name,method,matprop,TYPE,perimeter_case);
    [P_operator] = compute_P_operator(Msmooth,element,dim);
    regularization = 1;

    element.smoothing = 1;
    
    problem2solve.fhplots = fhplots;
end

element.Vfrac = Vfrac;
element.material.Vfrac = Vfrac;
element.material.Perimeter_target = Ptarget;
volume_integration = 'regularized_integration'; %regularized_integration naturalized_integration
ritz_gradient_representation_case = 'H1'; %'L2'; %only used for gradient method

global iter
iter = 0;

epsilon_scalar_product_P1 = 1*estimate_mesh_size(element,coordinates);
epsilon_iter_perimeter = create_epsilon_sequence(coordinates,element); 
epsilon_P0_gradient = 1*estimate_mesh_size(element,coordinates);
epsilon_Le_kernel = 1*estimate_mesh_size(element,coordinates);
problem2solve.scalar_product = @(f,g) f'*(((epsilon_scalar_product_P1)^2)*Stiff_smooth+Msmooth)*g;

[gamma0,gamma0_nodal,design_variable0] = ini_design_variables(dim,element,coordinates,Msmooth,problembsc,Stiff_smooth,P_operator,Algorithm,file_name,problem2solve.scalar_product,init_design_variable);
problem2solve.x0 = design_variable0;
problem2solve.gamma_hC0 = element.gamma_plus*ones(dim.nelem,1); %gamma0;
problem2solve.x0_elemental = gamma0;

[~,A_nodal_2_gauss,~,dvolu] = interpol(gamma0,element,dim,problembsc,coordinates);
M0 = sparse(1:dim.nelem,1:dim.nelem,dvolu);
% M0 = sparse(diag(dvolu));
smoothing_function = @(f_gp) smooth(dim.nelem,dim.npnod,dim.ndime,dim.nnode,f_gp,Msmooth,element,coordinates,coordinates,problembsc);

problem2solve.method = method;
problem2solve.algorithm = Algorithm;
problem2solve.element = element;
problem2solve.coordinates = coordinates;
problem2solve.problembsc = problembsc;
problem2solve.matprop = matprop;
problem2solve.problem_type = problem_type;
problem2solve.incropt = incropt;
problem2solve.design_type = TYPE;

problem2solve.P_operator = P_operator;
problem2solve.M0 = M0; 
problem2solve.Msmooth = Msmooth;
problem2solve.regularization = element.regularization;
problem2solve.epsilon_Le_kernel = epsilon_Le_kernel;

problem2solve.diff_react_equation = @(gamma,epsilon,kernel_case,variable_case) regularization_diff_reaction_equation(gamma,epsilon,kernel_case,Msmooth,Stiff_smooth,regularization,coordinates,element,problembsc,dim,P_operator,A_nodal_2_gauss,variable_case,M0,smoothing_function,problem2solve.algorithm);
ritz_gradient_representation = @(gradient) ritz_representation_of_gradient(gradient,epsilon_P0_gradient,Msmooth,Stiff_smooth,problem2solve.algorithm,ritz_gradient_representation_case);

problem2solve.epsilon_iter_perimeter =  epsilon_iter_perimeter;

% % % % % % % % % % % % % % % % % % % % %

% Optimization functions
problem2solve.equilibrium = @(gamma_gp,problembsc,h_C_0) equilibrium_update(gamma_gp,element,problembsc,fixnodes,coordinates,fext,fext_adjoint,dim,Msmooth,M0,Stiff_smooth,emass,h_C_0);
problem2solve.equilibrium_simp_all = @(gamma_gp,problembsc,element,h_C_0) equilibrium_update(gamma_gp,element,problembsc,fixnodes,coordinates,fext,fext_adjoint,dim,Msmooth,M0,Stiff_smooth,emass,h_C_0);

problem2solve.compliance = @(phi,problembsc,h_C_0) compliance_function(phi,epsilon_Le_kernel,kernel_case,problem2solve.equilibrium,problem2solve.diff_react_equation,ritz_gradient_representation,problembsc,Msmooth,h_C_0);
problem2solve.volume = @(phi,Vfrac) constraint_volume(phi,Vfrac,dim,element,problembsc,coordinates,Msmooth,M0,problem2solve.algorithm,volume_integration,epsilon_Le_kernel,kernel_case,problem2solve.diff_react_equation,ritz_gradient_representation);
problem2solve.perimeter = @(phi,epsilon) constraint_Perimeter(phi,epsilon,problem2solve.algorithm,Msmooth,M0,element.material.Perimeter_target,element,coordinates,problem2solve.diff_react_equation,ritz_gradient_representation);
problem2solve.enforceChdiff_square = @(phi,epsilon,compliance,typeCdifference,component,initial_values) enforce_Ch_difference_square_component(phi,epsilon,compliance,epsilon_Le_kernel,kernel_case,problem2solve.diff_react_equation,ritz_gradient_representation,Msmooth,typeCdifference,component,initial_values);
problem2solve.enforceChdiff = @(phi,compliance,typeCdifference,component,initial_values) enforce_Ch_difference_component(phi,compliance,epsilon_Le_kernel,kernel_case,problem2solve.diff_react_equation,ritz_gradient_representation,Msmooth,typeCdifference,component,initial_values);

% Numerical options
problem2solve.lb = zeros(size(gamma0_nodal));
problem2solve.ub = ones(size(gamma0_nodal));

% Print functions
print_function = @(iterplot,gamma_gp,gamma_nodal,gradient,phi,structural_values) printinfo(problembsc,file_write,iterplot,coordinates,element,dim,gamma_gp,gamma_nodal,gradient,phi,structural_values);
verbosity = 1; % 1/0
print_each_iter = 1; % 1/0
plot_each_iter = 1; % 1/0
problem2solve.outputopt = @(design_variable,optdata,last_iter) plot_data_opt(design_variable,optdata,print_function,TYPE,verbosity,print_each_iter,plot_each_iter,last_iter,select_screen);

% % % % % % % % % % % % % % % % % % % % % % % 
    
[xsol,fsol] = top_opt_algorithm(problem2solve);
                                         
end

