function [DOFr,dR] = DirichletCONDtime_homogZERO(DIRICHLET,DATA,ndim,MESH,GEOproperties,DATALOC)
% Goal. Determine DOFr and    dR(t) . Boundary conditions ZERO
% FLUCTUATIONS. HOMOGENIZATION
%
% JAHO- 20-MARCH-2O21
%--------------------------------------------------------------------------------------
if nargin == 0
    load('tmp.mat')
elseif nargin == 5
    DATALOC = [] ;
end

try
    [NODES_FACES,NODES_LINES] = NodesFacesLinesGID(DATA.NameFileMeshDATA(1:end-4)) ;
catch
    NAME_DATA = [NameFileMesh(1:end-4),'.gid',filesep,NameFileMesh(1:end-4)] ;
    [NODES_FACES,NODES_LINES] = NodesFacesLinesGID(NAME_DATA) ;
end

if size(MESH.COOR,2) ==2
    rnodALL = unique(cell2mat(NODES_LINES')) ;
else
    rnodALL = unique(cell2mat(NODES_FACES')) ;
end
MACRODEF = DIRICHLET.MACRODEF ;

nnode = length(rnodALL) ;

nloads = length(MACRODEF) ;
U =  zeros(nnode*ndim,2*nloads) ;  % there are two modes for each loading state 
a = zeros(2*nloads,length(DATA.STEPS)) ;
COORrel = MESH.COOR(rnodALL,:)' ;



for  iload = 1:nloads
    GRADuMACRO_fin = MACRODEF(iload).AMPLITUDE ;
    if iload == 1
        GRADuMACRO_ini = zeros(size(GRADuMACRO_fin)) ; 
    else
        GRADuMACRO_ini = MACRODEF(iload-1).AMPLITUDE ;
    end
    %  F_ident = Fgrad-eye(ndim) ;
    dMACROfin = GRADuMACRO_fin*COORrel ;
    dMACROfin = dMACROfin(:) ;  
    
    dMACROini = GRADuMACRO_ini*COORrel ;
    dMACROini= dMACROini(:) ;
    LIMITS_interval = MACRODEF(iload).INTERVAL ; 
      
    % Constant state (initial one )
    iload_ini = 2*iload-1 ;     
    U(:,iload_ini) =dMACROini ;    
    FactorSteps =  ones(size(DATA.STEPS)) ;
    FactorSteps = FactorSteps.*(DATA.STEPS >= LIMITS_interval(1)) ;
    if iload <nloads
    FactorSteps = FactorSteps.*(DATA.STEPS < LIMITS_interval(2)) ;
    else
        FactorSteps = FactorSteps.*(DATA.STEPS <= LIMITS_interval(2)) ;
    end
    a(iload_ini,:) = FactorSteps ; 
    % INCREMENT 
    iload_fin = 2*iload  ;     
    U(:,iload_fin) =dMACROfin-dMACROini ;    
    FactorSteps =  MACRODEF(iload).TIMEFUN_s(DATA.STEPS) ;
    FactorSteps = FactorSteps.*(DATA.STEPS >= LIMITS_interval(1)) ;
    if iload <nloads
    FactorSteps = FactorSteps.*(DATA.STEPS < LIMITS_interval(2)) ;
    else
        FactorSteps = FactorSteps.*(DATA.STEPS <= LIMITS_interval(2)) ;
    end
    a(iload_fin,:) = FactorSteps ;   
end

dR.U = U ;
dR.a = a ;

REPRESENT_PLOT_NORM =0 ; 

if REPRESENT_PLOT_NORM ==1 
    DISP = U*a; 
    ndisp = sqrt(sum(DISP.^2,1)) ;
    figure(919)
    hold on
    xlabel('Time step') 
    ylabel('Norm(dR)')
    plot(ndisp)
    
end





DOFr = small2large(rnodALL,ndim) ;  % Candidates for being DOFr





