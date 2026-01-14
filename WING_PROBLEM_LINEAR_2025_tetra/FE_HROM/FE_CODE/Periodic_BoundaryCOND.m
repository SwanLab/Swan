function  [DOFr DOFm G dR dispMACRO DATA] = Periodic_BoundaryCOND(COOR,CN,MACRODEF,CNb,DATA)
% Automatic enforcement of  Periodic Boundary conditions on   unit cells
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
% % Joaquín A. Hernández (jhortega@cimne.upc.edu), 29-Oct-2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%dbstop('21')
if nargin==0
    load('tmp1.mat')
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
% LINES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NODESln = DetermineLinesPeriodic_OLD(NODESpl) ;
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




%TOL_PERIODIC = 10 ;
%warning('Manually setting tolerance...')
[MASTER SLAVES] =MasterSlavesSets(NODESpl,NODESln,NODESpnt,COOR,TOL) ;

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
[Gi nods  nodm DIFFcoor]= GmatrixSlaveMaster(MASTER,SLAVES,COORref')  ;

%Assembly of G
G= sparse(3*size(Gi,1),3*size(Gi,2)) ;
n = size(Gi,1) ; m=size(Gi,2) ;
for i=1:3
    ROWS = i:3:3*n ; 
    COLS = i:3:3*m ; 
    G(ROWS,COLS ) = Gi ;
end

% Prescribed displacement
% % 3) Displacement are measured taking as reference NODESpnt(1)
% %  Accordingly, we redefine the coordinate matrix
% 2) Center of the unit cell


% Thus
uB = MACRODEF*DIFFcoor';

% dispMACRO
%dbstop('99')
dispMACRO = MACRODEF*COORref ;
dispMACRO = dispMACRO(:); 



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










% %
%
% rnod = {} ; uPRES={}  ;
% % 1) List of nodes at which displacement is prescribed  in DIRECTION  i = 1
% % All nodes pertaining to plane z = 0
% % Boundary nodes
% BoundaryNodes= unique(CNb(:)) ;
% xmin = min(COOR(BoundaryNodes,1)) ; xmin = xmin(1) ;
% idim=1 ;
% rnodBASEloc = find(abs(COOR(BoundaryNodes,1)-xmin)<1e-10) ;
% rnodBASE = BoundaryNodes(rnodBASEloc) ;
% rnod{idim} =rnodBASE;
% % Vector of prescribed displacements
% displ1 = 0 ;
% uPRES{idim} = displ1*ones(size(rnod{idim})) ;
% % 2) List of nodes at which displacement is prescribed  in DIRECTION  i = 2
% idim = 2;
% rnod{idim} =rnodBASE ;
% % Vector of prescribed displacements
% displ2 = 0 ;
% uPRES{idim} = displ2*ones(size(rnod{idim})) ;
% % 2) List of nodes at which displacement is prescribed  in DIRECTION  i = 3
% idim = 3;
% rnod{idim} =rnodBASE ;
% % Vector of prescribed displacements
% displ3 = 0 ;
% uPRES{idim} = displ3*ones(size(rnod{idim})) ;
% %%%% Set of restricted degrees of freedom and vector of prescribed
% %%%% displacements (dR)
% DOFr = [] ; dR = [] ;
% for idim = 1:ndim
%     DOFr = [DOFr ; (rnod{idim}-1)*ndim+idim];
%     dR = [dR ; uPRES{idim}];
% end