function  [DOFr DOFm G dR dispMACRO NODESpl] = Periodic_BoundaryCOND_lam(COOR,CN,strainINP,CNb,DATA)
% Automatic enforcement of  Periodic Boundary conditions on   unit cells
% -------------------------------------------------------------
% INPUT:
% --------------------------------------------
% COOR: Matrix of coordinates
% CN:  Connectivity matrix (volume)
% CNb: Connectivity matrix (boundary)
% strainINP: Macroscopic deformation
% ------------------------------------------------------------
% OUTPUT:
% DOFr: Slave DOFs
% DOFm: Master DOFs
% dR and G defined the relation between  d(DOFr) and d(DOFm)
%  d(DOFr) = dR + G*d(DOFm)
% dispMACRO: strainMACRO*COORref
% ------------------------------------------------------------------
% % Joaquín A. Hernández (jhortega@cimne.upc.edu), 26-Nov-2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%dbstop('21')
if nargin==0
    load('tmp1.mat')
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
h = zmax-zmin ; 
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
TOL = ChooseTolerance(CN,COOR) ;
% ---------------------------------
% POINTS PERTAINING TO  PLANES XMAX=0, XMIN= 0....
% --------------------
NODESpl =  DetermineePlanesPeriodic(COORbound,TOL,xmax,xmin,ymax,ymin,zmax,zmin,NODESbound) ;
% -------------------------------------------------
% LINES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NODESln = DetermineLinesPeriodic(NODESpl) ;
%%% POINTS
NODESpnt = DeterminePointsPeriodic(NODESpl) ;

%%% REMOVE INTERSECTIONS
% From planes
NODESpl = RemoveLinesPlanes(NODESpl,NODESln) ;
% From lines
NODESln = RemovePointsLine(NODESln,NODESpnt) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   MASTER/SLAVE NODES
% -----------------------------------------------------





[MASTER SLAVES] =MasterSlavesSets_lam(NODESpl,NODESln,NODESpnt,COOR,TOL) ;

%% In a laminate indNODMASTER{3}  is no longer  master

%%

%%%%%%%%%%%%%%%%%%%%%%%%%
% Prescribed displacement

% %  Accordingly, we redefine the coordinate matrix
% 2) Center of the unit cell
%

switch DATA.REFERENCE_POINT
    case 'CORNER'
        C = COOR(NODESpnt(1),:)';
        
    case 'CENTER'  % % % 3) Displacement are measured taking as reference NODESpnt(1)
        C = (DATA.CENTER)' ;
    otherwise
        error('OPTION NOT IMPLEMENTED')
end
COORref = (bsxfun(@plus,COOR',-C));
% Assembly of matrix Gi
%----------------------
[Gi nods  nodm COORsm COORcornerREF]= GmatrixSlaveMasterNEW(MASTER,SLAVES,COORref')  ;

%  
dB.SLAVE= CoarseDisplace(strainINP,COORsm.SLAVES',DATA,h,[]) ;
dB.MASTER= CoarseDisplace(strainINP,COORsm.MASTER',DATA,h,[]) ;
uB = dB.SLAVE-dB.MASTER ; 

 % Reference corner displacmeent 
dispREF =  CoarseDisplace(strainINP,COORcornerREF',DATA,h,[]) ;

dispMACRO= CoarseDisplace(strainINP,COORref,DATA,h,dispREF) ;





%dispMACRO = bsxfun(@plus,dispMACRO,-dispREF) ; 

% 
% Computation of displacements (without integration constant)

%%%% % Displacement of a given corner is set to zero 
% corner  = NODESpnt(1) ;
% dispCORNER = dispMACRO(:,corner) ; 
% % Therefore, the total displacement is 
% dispMACRO = bsxfun(@plus,dispMACRO,-dispCORNER) ; 
% Displacement of all points 

 dispMACRO = dispMACRO(:);
% 





%Assembly of G
G= sparse(3*size(Gi,1),3*size(Gi,2)) ;
n = size(Gi,1) ; m=size(Gi,2) ;
for i=1:3
    ROWS = i:3:3*n ;
    COLS = i:3:3*m ;
    G(ROWS,COLS ) = Gi ;
end
 

DOFr = zeros(3*length(nods),1) ;
dR = zeros(3*length(nods),1) ;
DOFm = zeros(3*length(nodm),1) ;
for i=1:3
    % dR = [dR; uB(i,:)'];
    ROWS = i:3:3*length(nods) ;
    DOFr(ROWS) = [  3*(nods-1)+i];
    dR(ROWS) = uB(i,:)' ;
    ROWS = i:3:3*length(nodm) ;
    DOFm(ROWS) = [ 3*(nodm-1)+i];
    
end

