function  [DOFr DOFm G dR dispMACRO DATA] = ZERO_KBCS_cub(COOR,CN,MACRODEF,CNb,DATA,TypeElementB)
% Automatic enforcement of  Zero Boundary conditions on   unit cells
% -------------------------------------------------------------
% INPUT:
% --------------------------------------------
% COOR: Matrix of coordinates
% CN:  Connectivity matrix (volume)
% CNb: Connectivity matrix (boundary)
% MACRODEF: Macroscopic deformation
% ------------------------------------------------------------
% OUTPUT:
% DOFr: Slave DOFs
% DOFm: Master DOFs
% dR and G defined the relation between  d(DOFr) and d(DOFm)
%  d(DOFr) = dR + G*d(DOFm)
% dispMACRO: strainMACRO*COORref
% ------------------------------------------------------------------
% % Joaquín A. Hernández (jhortega@cimne.upc.edu), 5-MAy-2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%dbstop('21')
if nargin==0
    load('tmp.mat')
    DATA.TOL_deter_BOUNDARYNODES  = 0.001; 
end

% 1) Determine ximin,ximax
NODESbound = unique(CNb(:)); % All boundary nodes
COORbound = COOR(NODESbound,:) ;
xmin = min(COORbound(:,1)) ; xmin = xmin(1) ;
xmax = max(COORbound(:,1)) ; xmax = xmax(1) ;
ymin = min(COORbound(:,2)) ; ymin = ymin(1) ;
ymax = max(COORbound(:,2)) ; ymax = ymax(1) ;
zmin = min(COORbound(:,3)) ; zmin = zmin(1) ;
zmax = max(COORbound(:,3)) ; zmax = zmax(1) ;

DATA.VOL_RVE = (xmax-xmin)*(ymax-ymin)*(zmax-zmin) ; 

xC = 0.5*(xmin+xmax); 
yC = 0.5*(ymin+ymax); 
zC = 0.5*(zmin+zmax); 
DATA.CENTER = [xC,yC,zC]; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4) Displacements induced by the macroscopic deformation (at boundary nodes)
%uMACRO = MACRODEF*COORbound';
% -------------------------------------------

% 5) Setting "master nodes" and "slaves nodes"
% -------------------------
% -----------------------------------------------
% Choosing tolerance
DATA = DefaultField(DATA,'TOL_deter_BOUNDARYNODES',[]) ; %.TOL_deter_BOUNDARYNODES
if isempty(DATA.TOL_deter_BOUNDARYNODES) 
TOL = ChooseTolerance(CN,COOR) ;
else 
    TOL = DATA.TOL_deter_BOUNDARYNODES ; 
end
%warning('Prescribing tolerance to ...')
%TOL = 10
% ---------------------------------
% POINTS PERTAINING TO  PLANES XMAX=0, XMIN= 0....
% --------------------
%TOL = 1e-10
NODESpl =  DetermineePlanesPeriodic(COORbound,TOL,xmax,xmin,ymax,ymin,zmax,zmin,NODESbound) ;
% -------------------------------------------------
% % LINES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NODESln = DetermineLinesPeriodic(NODESpl) ;
% %%% POINTS
% NODESpnt = DeterminePointsPeriodic(NODESpl) ;
% 

 
% Computation of matrix R  (nstrain x nelemBND*ndime)
NODESbnd = [] ;
for i = 1:size(NODESpl,1)
    for j = 1:size(NODESpl,2)
        NODESbnd = [NODESbnd;NODESpl{i,j}] ;
    end
end
NODESbnd = unique(NODESbnd) ;



%% MASTER/slaves DOFs  %
%dbstop('32')
ndim = size(COOR,2) ;
DOFr = Nod2DOF(NODESbnd,ndim) ; % Boundary DOFs
DOFm = [] ; 
G = [] ; 

  

  
% Prescribed displacement 
% -----------------------
switch DATA.REFERENCE_POINT
    case 'CORNER'
        C = COOR(NODESpnt(1),:)';
        
    case 'CENTER'  % % % 3) Displacement are measured taking as reference NODESpnt(1)
        C = (DATA.CENTER)' ;
    otherwise
        error('OPTION NOT IMPLEMENTED')
end

 
COORref =(bsxfun(@plus,COOR',-C));
%%%

% uBARp
% -----
dispMACRO = MACRODEF*COORref ;   % Macroscopic displacements of all DOFs
dispMACRO = dispMACRO(:); 
dR = dispMACRO(DOFr) ; 
% uBARs 
 



% %%% REMOVE INTERSECTIONS
% % From planes
% NODESpl = RemoveLinesPlanes(NODESpl,NODESln) ;
% % From lines
% NODESln = RemovePointsLine(NODESln,NODESpnt) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   MASTER/SLAVE NODES
% -----------------------------------------------------

% MASTER NODES
%%%%%%%%%%%%%%
 
% 
% 
% 
% %TOL_PERIODIC = 10 ;
% %warning('Manually setting tolerance...')
% [DOFr,DOFm,G,dR,dispMACRO] =     MasterSlavesMin(NODESpl,NODESpnt,CNb,COOR,TypeElementB,MACRODEF,DATA)  ;
% %%%%%%%%%%%%%%%%%%%%%%%%%
% % Prescribed displacement

% %  Accordingly, we redefine the coordinate matrix
% 2) Center of the unit cell
%