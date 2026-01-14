function  [DOFr DOFm G dR dispMACRO NODESpl] = ZERO_minTB_cal(COOR,CN,strainINP,CNb,DATA,TypeElementB)
% Automatic enforcement of  Periodic Boundary conditions + MINIMAL , on   laminate cells
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
%dbstop('23')
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
NODESzTOP = NODESpl{3,1} ; 
NODESzBOT = NODESpl{3,2} ; 

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
% Mastere/Slaves for periodic boundary conditions (all entities but zmin,zmax)
[ SLAVES ] =MasterSlavesSets_Zero(NODESpl,NODESln,NODESpnt,COOR,TOL) ;
% ----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Master/Slaves for minimal restrictions on top/bottom faces
disp('Calculating restrictions (mininal b.c.)')
[MASTER_uvFL,SLAVE_uvFL,G_uvFL,uBmacro] =...
    MasterSlavesMinTopBot(NODESpl,CNb,COOR,TypeElementB,strainINP,DATA,h,NODESzBOT,NODESzTOP) ; 
disp('Done')
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
%[Gi nods  nodm COORsm COORcornerREF]= GmatrixSlaveMasterNEW(MASTER,SLAVES,COORref')  ;


 % Reference corner displacmeent 
 corner = NODESpnt(1) ; 
 COORcornerREF = COOR(corner,:) ; 
 %
dispREF =  CoarseDisplace(strainINP,COORcornerREF',DATA,h,[]) ;
 COOR_slv = COOR(SLAVES,:) ; 
uB= CoarseDisplace(strainINP,COOR_slv',DATA,h,dispREF) ;
 nods = SLAVES ; 

%%%%%
%warning('Prueba')
%dB = dB_w ; uB = uB_w ;  nods = nods_w ; nodm = nodm_w ; Gi = Gi_w ; 
%  



%%% MACRO DISPLACEMENTS FOR ALL NODES 
% -------------------------------------------------------------
dispMACRO= CoarseDisplace(strainINP,COORref,DATA,h,dispREF) ;
dispMACRO = dispMACRO(:);
% ----------------------------------------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Assembly of vector of prescribed displacements and G
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[MASTER_uvFL,SLAVE_uvFL,G_uvFL] = MasterSlavesMinTopBot(NODESpl,CNb,COOR,TypeElementB);

% %------------------
nDOFr = 3*length(nods)   +2*length(SLAVE_uvFL);  % Number of slave DOFs 
nDOFm =  2*length(MASTER_uvFL)  ;  % Number of master DOFm

% -----------
% SLAVES DOFs 
% -----------
[DOFr dR ROWSslv] = AssemDOFslv_PERIOD_TBmin(nDOFr,nods,uB,SLAVE_uvFL,uBmacro) ; 
% -----------
% MASTER DOFs 
% -----------
[DOFm ROWSmst] = AssemDOFmst_Zero_TBmin(nDOFm,MASTER_uvFL) ; 
% Matrix G 
% --------
G= sparse(nDOFr,nDOFm) ;
%%% u, direction , periodic
%G(ROWSslv.u,ROWSmst.u) = Gi ; 
%%% v, direction  , periodic
%G(ROWSslv.v,ROWSmst.v) = Gi ; 
%%% w, direction  , periodic
%G(ROWSslv.w,ROWSmst.w) = Gi ; 
%%%% Additional constraints 
%--------------------------------
G(ROWSslv.uMIN,ROWSmst.uMIN) = G_uvFL ; 
%--------------------------------
G(ROWSslv.vMIN,ROWSmst.vMIN) = G_uvFL ; 

% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Now we sort DOFr and DOFm in ascending order
[DOFr iSLV] = sort(DOFr,'ascend') ; 
[DOFm iMST] = sort(DOFm,'ascend') ; 
% Therefore 
dR = dR(iSLV) ; 
G =G(iSLV,iMST) ; 

 


