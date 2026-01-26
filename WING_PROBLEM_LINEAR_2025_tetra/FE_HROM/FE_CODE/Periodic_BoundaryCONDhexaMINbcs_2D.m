function  [DOFr DOFm G dR dispMACRO DATA] = Periodic_BoundaryCONDhexaMINbcs_2D(COOR,CN,MACRODEF,CNb,DATA)
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
% See ALEX_THESIS_mine.pdf
% % Joaquín A. Hernández (jhortega@cimne.upc.edu), 29-Oct-2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%dbstop('22')
if nargin==0
    load('tmp0.mat')
    %  DATA.TOL_deter_BOUNDARYNODES  = 0.001;
end

% 1) Determine ximin,ximax
NODESbound = unique(CNb(:)); % All boundary nodes
COORbound = COOR(NODESbound,:) ;
xmin = min(COORbound(:,1)) ; xmin = xmin(1) ;
xmax = max(COORbound(:,1)) ; xmax = xmax(1) ;
ymin = min(COORbound(:,2)) ; ymin = ymin(1) ;
ymax = max(COORbound(:,2)) ; ymax = ymax(1) ;
zmax = [] ; zmin =[] ;

DATA.VOL_RVE = (xmax-xmin)*(ymax-ymin) ;

xC = 0.5*(xmin+xmax);
yC = 0.5*(ymin+ymax);

DATA.CENTER = [xC,yC];

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
% POINTS PERTAINING TO  LINES XMAX=0, XMIN= 0....
% --------------------
%TOL = 1e-10

[NODESln NormalPlanes]=  DetermineePlanesPeriodicHexaBCS(COORbound,TOL,xmax,xmin,ymax,ymin,zmax,zmin,NODESbound) ;
% -------------------------------------------------
%%% POINTS

lineamaster = 1 ;
lineaslave = [3] ;
NODOref = NODESln{lineamaster}(1) ;
NODESpnt=NODOref;
for islave=1:length(lineaslave)
    COORslaveREF = COOR(NODOref,:)' + NormalPlanes{lineamaster} ;
    COORslv = COOR(NODESln{lineaslave(islave)},:) ;
    % coorLOC = repmat(coorLOC,length(Nmast),1) ;
    diffCOOR = bsxfun(@minus,COORslv',COORslaveREF) ;
    normDIFF = sum(diffCOOR.*diffCOOR,1) ;
    [MINIM_c ,indiceF ]= min(normDIFF); % find(normDIFF<TOL) ;
    NODOrefSLV =    NODESln{lineaslave(islave)}(indiceF) ;
    
    NODESpnt = [NODESpnt; NODOrefSLV] ;
end







%%% REMOVE INTERSECTIONS
% From lines
NODESln = RemovePointsLinesHEXA(NODESln,NODESpnt) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   MASTER/SLAVE NODES
% -----------------------------------------------------




%TOL_PERIODIC = 10 ;
%warning('Manually setting tolerance...')
[MASTER SLAVES] =MasterSlavesSetsHexaBCS_2D(NODESln,NODESpnt,COOR,TOL,NormalPlanes) ;

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
ndim = size(COOR,2) ;
%Assembly of G
G= sparse(ndim*size(Gi,1),ndim*size(Gi,2)) ;
n = size(Gi,1) ; m=size(Gi,2) ;
for i=1:ndim
    ROWS = i:ndim:ndim*n ;
    COLS = i:ndim:ndim*m ;
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



DOFr = zeros(ndim*length(nods),1) ;
dR = zeros(ndim*length(nods),1) ;
DOFm = zeros(ndim*length(nodm),1) ;
for i=1:ndim
    % dR = [dR; uB(i,:)'];
    ROWS = i:ndim:ndim*length(nods) ;
    DOFr(ROWS) = [  ndim*(nods-1)+i];
    dR(ROWS) = uB(i,:)' ;
    ROWS = i:ndim:ndim*length(nodm) ;
    DOFm(ROWS) = [ ndim*(nodm-1)+i];
    
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