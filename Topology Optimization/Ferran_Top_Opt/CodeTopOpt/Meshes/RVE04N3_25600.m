% -----------------------------------------------
global problembsc
problembsc.problemtype = '2D';
problembsc.convergencecriteria = 'FORCE';
problembsc.convergenciavalue = 1e-6;
problembsc.linesearch_mecha = 'NO';
problembsc.totaltime = 1;   % segundos
problembsc.loadsteps = 1;
problembsc.deltat =problembsc.totaltime/problembsc.loadsteps;
problembsc.iterationsteps = 40;
problembsc.frqpostpros = 1;
problembsc.smoothing_proc = 2; % 1: lumped, 2:full matrix (afecta suavizado y calculo de theta)
problembsc.max_kappa_iter = 50;
problembsc.smoothDtC = 1;
problembsc.costfunct = 'MIN_INV_STIFF';%' 'MIN_MINUS_STIFF MIN_INV_STIFF MIN_STIFF_INV % HORIZONTAL,  BULK_MAX, SHEAR_MAX MAXSTIFF 
problembsc.restart = 'NO';
problembsc.meth_kappa_opt = 'COST_ONLY'; %IMPOSED BRUTE_FORCE FREQUENCY COST_MASS_THETA COST_ONLY COST_OPT
problembsc.lagrange_update = 'LINEAL' ; % AUGMENTED QUADRATIC

gravity = 9.810;
dyna.Gx = 0;       % debe estar mm/s2; se supone que las coord son 'mm'
dyna.Gy = -gravity;   % debe estar mm/s2;    


% Elements
%
global element
element.name = 'bloque';
element.type = 'TRIANGLE'; % TRIANGLE QUAD HEXA
element.formulation = 'SSSD_LT'; % LSLD_LT LSLD_LA
element.material.name = 'STEEL';
element.material.type = 'HOOKE_LAW'; 
element.material.subtype = 'PLANESTRES'; % PLANESTRAIN PLANESTRES
element.material.poiss = 0.3;
element.material.young = 1.0;
element.material.kappa  = element.material.young/(3*(1-2*element.material.poiss));
element.material.G      = element.material.young/(2*(1+element.material.poiss));
element.material.mu     = element.material.G;
element.material.lambda = element.material.young*element.material.poiss/...
                          ((1+element.material.poiss )*(1-2*element.material.poiss));
element.material.opt_epsi = 1.0e-2;
element.material.kappa_min = 1e-04; %bulk 1.0e-4;
element.material.kappa_ini = 0; % [ini,end] rango de estudio de  kappa
element.material.kappa_end = [1]; %bulk [0.5 0.15 0.05 0.005];
element.material.kappa_maxiter = 100;
element.material.kappa_reduction = 0.9; % 80%
element.material.vol_reduction = [20]; %bulk [15 5 2 2 1 1 1 0.5]; % en porcentaje
element.material.opt_L = 40; % bulk 20;
element.material.opt_eta = 0.0000000625;
element.material.thetaend = 1e-3; %theta mas peque�o que se pide, en grados
element.material.opt_regul = 1; % 0: ptos de gauss, 1: todo elemental
element.material.Emeth = 2; % 1: en gauss, 2: de E nodal, 3:cte por elemento (ver cal_young_modulus)   
element.material.Vfrac = 0.5;
element.material.penalty = 15;
% HISTORY OF KAPPA_OPT
% history of kappa used to obtain rve-horizontal with the coarse mesh (6400
% elem)
hkappa = 0.5*ones(1,problembsc.iterationsteps);
          

% history of kappa used to obtain rve-bulk with the coarse mesh (6400 elem)
%hkappa = [0.45 0.12 0.025 0.00253 1.22e-03 0.02];
% hkappa = [4.6e-001 4.6e-001 3.8e-002 2.4e-002 2.3e-002 1.97e-002 1.0e-002 1e-02 1.0e-002 1.0e-002...
%         2e-02 2e-02 2e-02 2e-02 2e-02 2e-02 2e-02 2e-02 2e-03 2e-03  ...
%         2e-03 2e-03 2e-03 2e-03 2e-03 2e-03 2e-03 2e-03 2e-03 2e-03  ... 
%         2e-03 2e-03 2e-03 2e-03 2e-03 2e-03 2e-03 2e-03 2e-04 2e-03  ...
%         2e-03 2e-03 2e-03 2e-03 2e-03 2e-03 2e-03 2e-03 2e-03 2e-03  ... 
%         2e-03 2e-03 2e-03 2e-03 2e-03 2e-03 2e-03 2e-03 2e-04 2e-03  ...
%         2e-03 2e-03 2e-03 2e-03 2e-03 2e-03 2e-03 2e-03 2e-03 2e-03  ... 
%         2e-03 2e-03 2e-03 2e-03 2e-03 2e-03 2e-03 2e-03 2e-04 2e-03  ...
%         2e-03 2e-03 2e-03 2e-03 2e-03 2e-03 2e-03 2e-03 2e-03 2e-03  ... 
%         2e-03 2e-03 2e-03 2e-03 2e-03 2e-03 2e-03 2e-03 2e-04 2e-03  ...
%         ];
element.material.hkappa = hkappa;    

% FICHERO CON LAS COORDENADAS Y CONECTIVIDADES
fname = 'RVE04N3_25600_MALLA'; 
eval(fname); 

nel=size(gidlnods,1);
ncol=size(gidlnods,2);
element.conectivities = gidlnods(1:nel,2:ncol);
npt=size(gidcoord,1);
ncol=size(gidcoord,2);
coordinates = gidcoord(1:npt,2:ncol-1); 
clear gidcoord;
clear gidlnods;

% Fixed Nodes
global fixnodes

%
fext.pointload = [] ;
%
% Side loads
%
fext.sideload = [ ];
% 


element.pnods = [];
element.lglib = [];
fixnodes = [];
BC = 'PERIODIC';
switch BC
    case 'PERIODIC'
        fname = 'RVE04N3_25600_PERIODIC';
        eval(fname);
        % prescription of the corners
        ifix=0;
        for i=1:size(corners,2)
            ifix=ifix+1;
            fixnodes(ifix,1)=corners(i); % node
            fixnodes(ifix,2)=1; % idim
            fixnodes(ifix,3)=0; % U_imp
            ifix=ifix+1;
            fixnodes(ifix,1)=corners(i); % node
            fixnodes(ifix,2)=2; % idim
            fixnodes(ifix,3)=0; % U_imp
        end

    case 'SIMPLE_APOYO'
        ifix=1; i = 3278;
        fixnodes(ifix,1)=i; % node
        fixnodes(ifix,2)=1; % idim
        fixnodes(ifix,3)=0; % U_imp
        ifix=ifix+1;
        fixnodes(ifix,1)=i; % node
        fixnodes(ifix,2)=2; % idim
        fixnodes(ifix,3)=0; % U_imp
        ifix=ifix+1; i = 4225;
        fixnodes(ifix,1)=i; % node
        fixnodes(ifix,2)=1; % idim
        fixnodes(ifix,3)=0; % U_imp
         ifix=ifix+1; i = 4225;
        fixnodes(ifix,1)=i; % node
        fixnodes(ifix,2)=2; % idim
        fixnodes(ifix,3)=0; % U_imp
    case 'TODO_IMPUESTO'
         ifix=0;fixnodes=[];
        h=0.015625;
        %vertical izq
        for i=1:size(coordinates,1)
            if (coordinates(i,1)<h/5)
                ifix=ifix+1;
                fixnodes(ifix,1)=i; % node
                fixnodes(ifix,2)=1; % idim
                fixnodes(ifix,3)=0; % U_imp
                ifix=ifix+1;
                fixnodes(ifix,1)=i; % node
                fixnodes(ifix,2)=2; % idim
                fixnodes(ifix,3)=0; % U_imp
            end
        end
        %vertical der
        for i=1:size(coordinates,1)
            if (coordinates(i,1)>1-h/5)
                ifix=ifix+1;
                fixnodes(ifix,1)=i; % node
                fixnodes(ifix,2)=1; % idim
                fixnodes(ifix,3)=0; % U_imp
                ifix=ifix+1;
                fixnodes(ifix,1)=i; % node
                fixnodes(ifix,2)=2; % idim
                fixnodes(ifix,3)=0; % U_imp
            end
        end
        %horizontal inferior
        for i=1:size(coordinates,1)
            if (coordinates(i,2)<h/5)
                ifix=ifix+1;
                fixnodes(ifix,1)=i; % node
                fixnodes(ifix,2)=1; % idim
                fixnodes(ifix,3)=0; % U_imp
                ifix=ifix+1;
                fixnodes(ifix,1)=i; % node
                fixnodes(ifix,2)=2; % idim
                fixnodes(ifix,3)=0; % U_imp
            end
        end
        %horizontal superior
        for i=1:size(coordinates,1)
            if (coordinates(i,2)>1-h/5)
                ifix=ifix+1;
                fixnodes(ifix,1)=i; % node
                fixnodes(ifix,2)=1; % idim
                fixnodes(ifix,3)=0; % U_imp
                ifix=ifix+1;
                fixnodes(ifix,1)=i; % node
                fixnodes(ifix,2)=2; % idim
                fixnodes(ifix,3)=0; % U_imp
            end
        end
end
element.nodesolid = [];

% INITIAL VALUE OF THE LEVEL SET FUNCTION
x = coordinates(:,1); x0 = 0.5*ones(size(coordinates,1),1);
y = coordinates(:,2); y0 = 0.5*ones(size(coordinates,1),1);
c2x = cos(pi*(x-x0)).^2;
c2y = cos(pi*(y-y0)).^2;
% N = sqrt(0.140474562963467), N del articulo que calculo numericamente
phifunct_n = 1/sqrt(0.140474562963467)*(c2x.*c2y - 0.5);

vol_n = 0.6;

%[ phifunct ] = ini_phifunct(size(coordinates,1));

problembsc.ppinfo(1) = 0; %CAUCHY GAUSS POINT
problembsc.ppinfo(2) = 0; % b_e gauus point
problembsc.ppinfo(3) = 0; % alpha 
problembsc.ppinfo(4) = 0; % J 
problembsc.ppinfo(5) = 0; % Von Mises gauss point
problembsc.ppinfo(6) = 0; 
problembsc.ppinfo(7) = 0; 
problembsc.ppinfo(8) = 0; % los desplazamientos de los tres experimentos   
problembsc.ppinfo(9) = 0; % tdisp 
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
problembsc.ppinfo(28) = 1; % despl total - mixed model
problembsc.ppinfo(29) = 0; % funcion g (pto de gauss,cuando no se usa gnodal)
problembsc.ppinfo(30) = 1; % funcion E modulus
problembsc.ppinfo(31) = 1; % funcion phi 
problembsc.ppinfo(32) = 1; % plot cost function and theta
problembsc.ppinfo(33) = 0; % funcion g_til (gauss a partir de gnodal)
problembsc.ppinfo(34) = 1; % funcion g orthogonal


