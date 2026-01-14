function  [DOFr DOFm G dR dispMACRO NODESpl] = Periodic_BoundaryCOND_lamMIN(COOR,CN,strainINP,CNb,DATA)
% Automatic enforcement of  Periodic Boundary conditions on   unit cells
% Version that ensures consistency with transverse shear strains 
% See document Plates_comp_homogenization.pdf 
%              ---------------------------------
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
% % Joaquín A. Hernández (jhortega@cimne.upc.edu), 18-Dic-2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%dbstop('24')
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
%NODESln = DetermineLinesPeriodic(NODESpl) ;
NODESln = DetermineLinesPeriodicOLD(NODESpl) ; % Modified 13-Dec-2020, copied from  
%/home/joaquin/Desktop/CURRENT_TASKS/COMPOSITE_MATERIALS_DOCENCIA/APUNTES_LATEX/FINITE_ELEMENT_CODE/FE_CODE/Periodic_BoundaryCOND_lamTB.m
 
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
% For disp. u,v --> Include all surfaces
[MASTER SLAVES] =MasterSlavesSets_lamTB(NODESpl,NODESln,NODESpnt,COOR,TOL) ;
% ----------------------------------------------------------------------------
% For disp. w --> Exclude set "i=3"
iZMAX = 3 ; 
nn = size(MASTER,2) ; mmm = size(SLAVES,2) ;
MASTER_w = cell(1,nn-1);
SLAVES_w = cell(nn-1,mmm) ; 
iii = 1 ; 
for i = 1:nn
    if i ~= iZMAX
        MASTER_w{iii} = MASTER{i} ; 
        for jjj=1:mmm
        SLAVES_w{iii,jjj} = SLAVES{i,jjj} ; 
        end
        iii = iii +1 ;             
    end
end
 
% ----------------------------------------------------------------------------
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
% Assembly of matrix Gi for u and v
%-----------------------------------
[Gi nods  nodm COORsm COORcornerREF]= GmatrixSlaveMasterNEW(MASTER,SLAVES,COORref')  ;
dB.SLAVE= CoarseDisplace(strainINP,COORsm.SLAVES',DATA,h,[]) ;
dB.MASTER= CoarseDisplace(strainINP,COORsm.MASTER',DATA,h,[]) ;
uB = dB.SLAVE-dB.MASTER ; 
% Assembly of matrix Gi for w
%-----------------------------------
[Gi_w nods_w  nodm_w COORsm_w COORcornerREF_w]= GmatrixSlaveMasterNEW(MASTER_w,SLAVES_w,COORref')  ;
dB_w.SLAVE= CoarseDisplace(strainINP,COORsm_w.SLAVES',DATA,h,[]) ;
dB_w.MASTER= CoarseDisplace(strainINP,COORsm_w.MASTER',DATA,h,[]) ;
uB_w = dB_w.SLAVE-dB_w.MASTER ; 

%%%%%
%warning('Prueba')
%dB = dB_w ; uB = uB_w ;  nods = nods_w ; nodm = nodm_w ; Gi = Gi_w ; 
%  


 % Reference corner displacmeent 
dispREF =  CoarseDisplace(strainINP,COORcornerREF',DATA,h,[]) ;
%%% MACRO DISPLACEMENTS FOR ALL NODES 
% -------------------------------------------------------------
dispMACRO= CoarseDisplace(strainINP,COORref,DATA,h,dispREF) ;
dispMACRO = dispMACRO(:);
% ----------------------------------------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Assembly of vector of prescribed displacements and G
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %------------------
nDOFr = 2*length(nods) + length(nods_w) ;  % Number of slave DOFs 
nDOFm = 2*length(nodm) + length(nodm_w) ;  % Number of master DOFm

% -----------
% SLAVES DOFs 
% -----------
[DOFr dR ROWSslv] = AssemDOFslv_PERIOD_TB(nDOFr,nods,nods_w,uB,uB_w) ; 
% -----------
% MASTER DOFs 
% -----------
[DOFm ROWSmst] = AssemDOFmst_PERIOD_TB(nDOFm,nodm,nodm_w) ; 
% Matrix G 
% --------
G= sparse(nDOFr,nDOFm) ;
%%% u, direction 
G(ROWSslv.u,ROWSmst.u) = Gi ; 
%%% v, direction 
G(ROWSslv.v,ROWSmst.v) = Gi ; 
%%% w, direction 
G(ROWSslv.w,ROWSmst.w) = Gi_w ; 
% 
%%%% Now we sort DOFr and DOFm in ascending order
[DOFr iSLV] = sort(DOFr,'ascend') ; 
[DOFm iMST] = sort(DOFm,'ascend') ; 
% Therefore 
dR = dR(iSLV) ; 
G =G(iSLV,iMST) ; 

 


