settings.ptype = 'MACRO';
% settings.filename = 'CantileverBeam_Triangle_Linear_Fine';
% settings.filename = 'topopt_quad';
% settings.filename = 'GrippingNew';
% settings.ptype = 'MICRO';
% settings.filename = 'RVE_Square_Triangle';
% settings.filename = 'RVE_Square_Triangle_Fine';

settings.method = 'SIMPALL';
%settings.method = 'SIMP_P3';
%settings.method = 'SIMP_Adaptative';

settings.material = 'ISOTROPIC';
settings.initial_case = 'full';
% settings.initial_case = 'circle';
% settings.initial_case = 'horizontal';
% settings.initial_case = 'square';
% settings.initial_case = 'feasible';
% settings.initial_case = 'rand';



settings.cost = {'compliance';'perimeter'}; %'chomog_fraction';'compliance';'perimeter';'chomog_alphabeta';'nonadjoint_compliance';
% settings.weights = []; %all 1
settings.weights = [1 0.1]; %compl+lambda*perimeter
settings.constraint = {'volumeConstraint'};


settings.optimizer = 'SLERP';
%settings.optimizer = 'PROJECTED GRADIENT';settings.incrementFactor = 1;
%settings.optimizer = 'MMA';
%settings.optimizer = 'IPOPT';


settings.filterType = 'P1';
%settings.filter = 'PDE';

settings.TOL.rho_plus = 1;
settings.TOL.rho_minus = 0;
settings.TOL.E_plus = 1;
settings.TOL.E_minus = 1e-3;
settings.TOL.nu_plus = 1/3;
settings.TOL.nu_minus = 1/3;



settings.nsteps = 1;
settings.Vfrac_final = settings.target_parameters.Vfrac;
settings.optimality_final = settings.target_parameters.optimality_tol;
settings.constr_final = settings.target_parameters.constr_tol;
settings.Vfrac_initial = 1;


settings.optimality_initial = 1e-3;
settings.constr_initial = 1e-3;


settings.micro.alpha = [1 1 0]';
settings.micro.beta = [1 1 0]';