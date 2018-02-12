function [fext,fext_adjoint,element,fixnodes,problembsc,coordinates,dim,Msmooth,Stiff_smooth,emass,Group] = read_data_problem(file_name,TYPE,perimeter_case)

% -----------------------------------------------
problembsc.problemtype = '2D';
problembsc.convergencecriteria = 'FORCE';
problembsc.convergenciavalue = 1e-6;
problembsc.linesearch_mecha = 'NO';
problembsc.totaltime = 1;   % segundos
problembsc.loadsteps = 1;
problembsc.deltat =problembsc.totaltime/problembsc.loadsteps;
problembsc.iterationsteps = 2000;%2000;%2000;%2000; %40
problembsc.frqpostpros = 1;
problembsc.smoothing_proc = 2; % 1: lumped, 2:full matrix (afecta suavizado y calculo de theta)
problembsc.max_kappa_iter = 50;
problembsc.smoothDtC = 1;
problembsc.smoothChi = 2; % 0: o 0 o 1 en txi gauss a traves de phi interpolado (gauss), 1: 0 1/3 2/3 1 interpolado por txi nodal, 2: continua, interp por phi 3:cte por elemento  
problembsc.costfunct = 'MIN_INV_STIFF';%' 'MIN_MINUS_STIFF MIN_INV_STIFF MIN_STIFF_INV % HORIZONTAL,  BULK_MAX, SHEAR_MAX MAXSTIFF 
problembsc.restart = 'NO';
problembsc.meth_kappa_opt = 'COST_OPT';%'COST_OPT';%'COST_BEST';%'COST_BEST';%'COST_OPT' 'COST_ONLY_INC'%'COST_ONLY_INC'%IMPOSED BRUTE_FORCE FREQUENCY COST_MASS_THETA COST_ONLY COST_OPT 'COST_ONLY_REMESH'
problembsc.lagrange_update = 'POTENCIAL' ; % AUGMENTED POTENCIAL
problembsc.constraint = 'consistent'; % no_consistent
problembsc.potential_constraint = 1;
problembsc.nkappa = 50;
problembsc.penalty_factor = 0.5;%0.99;%0.999;
problembsc.flag_change_micro = 0; % 1 si en el proceso macro, 0 no
problembsc.flag_change_macro = 1;
problembsc.TOL = 1e-3;
problembsc.macro_end_flag_by_groups = 0; %0 es punto de gauss, %1 es por grupos
problembsc.macro_on_flag_by_groups = 0; %0 es punto de gauss, %1 es por grupos
problembsc.macro_end_flag_postproces_micro = 1; % 1 si se postprocesa, 0 no
problembsc.macro_on_flag_postproces_micro = 0; % 1 si se postprocesa, 0 no
problembsc.mirroring = 1; %Post-proceso, if non symmetric, do mirroring (macroscale)
problembsc.algorithm_update = 'AMSTUTZ';%'BCN_phig_EX';%'BCN_phig_SI';%'BCN_gphi_EX';%'AMSTUTZ';%BCN_gphi_SI';%'BCN'; %'AMSTUTZ' 
problembsc.regularization_perimeter = 'YES';%'NO'
problembsc.augmented_type = 'NOVOTNY';%'NOVOTNY';%'BIS'; %'NOCEDAL'; 'THETA_DEPENDANCY'; %'INCREASE','BIN'
problembsc.phisical_type = 'ELASTIC';%ELASTIC';'THERMAL';%ELASTIC';'THERMAL';%ELASTIC';'THERMAL';%ELASTIC';
problembsc.force_vol = 1;
problembsc.alpha_perimeter = 0;
problembsc.tol_epsi = 1e-6;
problembsc.max_epsilon_iter = 1;
problembsc.penalty_max = 10;
gravity = 9.810;
dyna.Gx = 0;       % debe estar mm/s2; se supone que las coord son 'mm'
dyna.Gy = -gravity;   % debe estar mm/s2;  


element.random_first_micro = 0; % O first micro homogeneos, 1 first micro randomnly heterogeneos
element.name = 'bloque';
element.type = 'TRIANGLE'; % TRIANGLE QUAD HEXA
element.material.name = 'STEEL';

switch problembsc.phisical_type
    case 'THERMAL'
        element.material.type = 'FOURIER_LAW'; 
    case 'ELASTIC'
        element.material.type = 'HOOKE_LAW'; 
end
element.material.subtype = 'PLANESTRES'; % PLANESTRAIN PLANESTRES
element.material.homogenous_material = 'YES';%'VADEMECUM'; %'YES' 'VADEMECUM'
element.material.anisotropic = 'ANISOTROPIC'; %'ISOTROPIC
element.material.poiss = 0.3;
element.material.young = 1.0;
element.material.k11 = 1;
element.material.k12 = 0;
element.material.k22 = 1;
element.material.kappa  = element.material.young/(3*(1-2*element.material.poiss));
element.material.G      = element.material.young/(2*(1+element.material.poiss));
element.material.mu     = element.material.G;

element.gamma_minus = 0;
element.gamma_plus = 1;

element.E_plus = 1;
element.E_minus = 0.001;
element.nu_plus = -0.6;
element.nu_minus = 0.33333;

Vfrac = 0.3;
Perimeter_target = 1;
element.Vfrac = Vfrac;


lambda = element.material.young*element.material.poiss/...
    ((1+element.material.poiss )*(1-2*element.material.poiss));
switch element.material.subtype
    
    case 'PLANESTRAIN'
        element.material.lambda = lambda;
    case  'PLANESTRES'
        mu = element.material.mu;
        element.material.lambda = 2*mu*lambda/(lambda+2*mu);
end

element.material.opt_epsi = 1.0e-6;
element.material.kappa_min = 0.005; %bulk 1.0e-4;
element.material.kappa_ini = 0; % [ini,end] rango de estudio de  kappa
element.material.kappa_end = [1]; %bulk [0.5 0.15 0.05 0.005];
element.material.kappa_maxiter = 100;
element.material.kappa_reduction = 0.5; % 80%
element.material.vol_reduction = [20]; %bulk [15 5 2 2 1 1 1 0.5]; % en porcentaje
element.material.opt_L = 40; % bulk 20; %100
element.material.opt_eta = 0.0000000625;
element.material.thetaend = 1; %theta mas peque�o que se pide, en grados
element.material.theta0 = 5;
element.material.opt_regul = 0.1; % 0: ptos de gauss, 1: todo elemental
element.material.Emeth = 0; % 0: en gauss, 1: de E nodal, 2:cte por elemento (ver cal_young_modulus)   
element.material.Vfrac = Vfrac;
element.material.Perimeter_target = Perimeter_target;
element.material.penalty = 10;
element.Vol_inc_tol = 5;%0.05;
element.fobj_inc_tol = 0;%0.1;%5;%1e2;
element.ngaus = 1;%3
element.regularization = 1;

Vol_ini = 0.7987;

hkappa = 0.5*ones(1,problembsc.iterationsteps);


element.material.hkappa = hkappa;    
pointload = [];
sideload = [];
nodesolid = [];
lnodes = [];
corners = [];
Group = [];
Initial_holes = [];
External_border_nodes = [];

% FICHERO CON LAS COORDENADAS Y CONECTIVIDADES
imesh = 1;
fname = [file_name,'_MESH_',num2str(imesh)];%'RVE04N3_MALLA'; 
eval(fname); 

if exist('Micro_slave','var')
    element.Micro_slave = Micro_slave;
end
    
nel=size(gidlnods,1);
if size(gidlnods,2) == 5   
       if gidlnods(1,5) == 0 ;
        ncol = size(gidlnods,2)-1;
        gidlnods = gidlnods(:,1:ncol);
       end
else
ncol=size(gidlnods,2);
end
element.conectivities = gidlnods(1:nel,2:ncol);
npt=size(gidcoord,1);
ncol=size(gidcoord,2);
coordinates = gidcoord(1:npt,2:ncol-1); 
clear gidcoord;
clear gidlnods;

% Pointload
if exist('pointload_complete','var')
    fext.pointload = pointload_complete ;
else
    fext.pointload = pointload ;
end

% Side loads
fext.sideload = sideload;

%Adjoint
if exist('pointload_adjoint','var')
fext_adjoint.pointload = pointload_adjoint;
else
fext_adjoint.pointload = [];    
end
%micro of gauss point to post process
element.micro_2_post = Micro_gauss_post;

%holes
if isempty(Initial_holes)
element.initial_holes = Initial_holes;
else
element.initial_holes = Initial_holes(:,1);
end

switch perimeter_case
    case 'domain_perimeter'
        element.External_border_nodes = [];
    case 'domain_and_contour_perimeter'
        if isempty(External_border_nodes) && strcmp(TYPE,'MACRO')
            error('This option is not available with this mesh');
        end
        element.External_border_nodes = [External_border_nodes , ones(size(External_border_nodes))];
end

% Initialize basics dimensions
[dim.npnod,dim.nndof,dim.ndime,dim.nunkn,dim.nstre] = data_nod(coordinates, element.type,element,problembsc);
[dim.nelem,dim.nnode,dim.neleq] = data_elem(element.conectivities, element.type);

[~,~,ngaus] = cal_posgp_weigp(element.type,dim.ndime,dim.nnode,element.ngaus);
element.ngaus = ngaus;
element.material.young = element.material.young;%*ones(1,dim.nelem);%rand(1,dim.nelem);

if isempty(pointload)
    boundary_nodes = [];
else
boundary_nodes = unique([pointload(:,1);lnodes(:,1)]);
end

element.boundary_elements = [];

element.dist_frontera = max(sqrt(coordinates(:,1).^2 + coordinates(:,2).^2))/100;
for ielem = 1:dim.nelem
    nodes_elem = element.conectivities(ielem,:);
    xbar = sum(coordinates(nodes_elem,1))/dim.nnode;
    ybar = sum(coordinates(nodes_elem,2))/dim.nnode;
    element.baricenter(ielem,:) = [xbar ybar];
end

element.boundary_elements = Boundary_elements;

element.boundary_nodes = boundary_nodes;
element.lglib = [];
problembsc.TYPE = TYPE;


switch problembsc.TYPE
    case 'MACRO'
        BC = 'NORMAL_DIRICLET';
        BC_perimeter = 'NULL_DIRICLET';
        [fixnodes,pnods] = compute_fix_nodes(BC,coordinates,lnodes,problembsc.phisical_type,file_name);
        [fixnodes_perimeter,pnods_perimeter] = compute_fix_nodes(BC_perimeter,coordinates,element.External_border_nodes,'THERMAL',file_name);
        nunkn_per = 1;
        [fix_df_per,free_df_per] = fix_free_degree_freedom(element.type,dim.npnod*nunkn_per,nunkn_per,fixnodes_perimeter,problembsc);
        element.fix_df_per = fix_df_per;
        element.free_df_per = free_df_per;
    
    case 'MICRO'
        BC = 'PERIODIC';
        BC_perimeter = 'PERIODIC';
        [fixnodes,pnods] = compute_fix_nodes(BC,coordinates,lnodes,problembsc.phisical_type,file_name);
        [fixnodes_perimeter,pnods_perimeter] = compute_fix_nodes(BC_perimeter,coordinates,element.External_border_nodes,'THERMAL',file_name);
        element.fix_df_per = [];
        element.free_df_per = [];
        element.fixnodes_perimeter = [];

end

%element.nodesolid = nodesolid;
element.pnods = pnods;
element.fixnodes = fixnodes;
element.pnods_perimeter = pnods_perimeter;
element.fixnodes_perimeter = fixnodes_perimeter;


%element.boundary_nodes = coordinates(:,1) == max(coordinates(:,1)) | coordinates(:,1) == min(coordinates(:,1)) | coordinates(:,2) == max(coordinates(:,2)) | coordinates(:,2) == min(coordinates(:,2));

% basic variables
[coordinatesn,coordinatesa] = init_coord(coordinates);
[Msmooth,emass] = mass_matrix(dim,problembsc,element,coordinatesn,coordinatesa);
[Stiff_smooth] = stiff_unitary_triang( dim,element,problembsc,coordinatesn,coordinatesa,2,1,dim.npnod*1);
%[ phifunct ] = ini_phifunct(size(coordinates,1));

if problembsc.force_vol 
[force_ext_vol] = fext_vol(dim,coordinatesn,coordinatesa,element,problembsc,coordinates);
end
fext.fvol = force_ext_vol;
fext_adjoint.fvol = 0*force_ext_vol;

problembsc.ppinfo(1) = 1; %CAUCHY GAUSS POINT
problembsc.ppinfo(2) = 0; % b_e gauus point
problembsc.ppinfo(3) = 0; % alpha 
problembsc.ppinfo(4) = 1; % J 
problembsc.ppinfo(5) = 0; % Von Mises gauss point
problembsc.ppinfo(6) = 0; 
problembsc.ppinfo(7) = 0; 
problembsc.ppinfo(8) = 0; % los desplazamientos de los tres experimentos   
problembsc.ppinfo(9) = 1; % tdisp 
problembsc.ppinfo(10) = 0; % velocity 
problembsc.ppinfo(11) = 0; % aceleration
problembsc.ppinfo(12) = 0; % nodal pressure - mixed model
problembsc.ppinfo(13) = 0; % inc_disp - mixed model
problembsc.ppinfo(14) = 0; % elastico o plastico - code
problembsc.ppinfo(15) = 0; % media kirchhoff stress
problembsc.ppinfo(16) = 0; % plastic disipation
problembsc.ppinfo(17) = 0; % STRAIN RATE
problembsc.ppinfo(18) = 0; % alpha factor = 1/(.. + lambda)
problembsc.ppinfo(19) = 0; % b1
problembsc.ppinfo(20) = 0; % Radios circunscrito e inscrito
problembsc.ppinfo(21) = 0; % fuerza residual
problembsc.ppinfo(22) = 0; % fuerza inercial
problembsc.ppinfo(23) = 0; % fuerza interna
problembsc.ppinfo(24) = 0; % fuerza interna_n
problembsc.ppinfo(25) = 0; % fuerza volumen
problembsc.ppinfo(26) = 0; % media CAUCHY stress
problembsc.ppinfo(27) = 0; % fuerza de contacto
problembsc.ppinfo(28) = 0; % despl total - mixed model
problembsc.ppinfo(29) = 0; % funcion g (pto de gauss,cuando no se usa gnodal)
problembsc.ppinfo(30) = 1; % funcion E modulus
problembsc.ppinfo(31) = 1; % funcion phi 
problembsc.ppinfo(32) = 1; % plot cost function and theta
problembsc.ppinfo(33) = 1; % funcion g_til (gauss a partir de gnodal)
problembsc.ppinfo(34) = 1; % funcion g orthogonal
problembsc.ppinfo(35) = 0; % funcion phi y theta de la micro
problembsc.ppinfo(36) = 1; % Ch
problembsc.ppinfo(37) = 0; % Fobj local
problembsc.ppinfo(38) = 0; % Group
end


