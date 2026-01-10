function  [DISP_CONDITIONS,COORref,DATA]  = Periodic_BoundaryCOND_LARGE(COOR,CN,CNb,DATA)
%==========================================================================
% Periodic_BoundaryCOND_LARGE
%
% PURPOSE
%   Automatic construction of **periodic boundary conditions (PBCs)** for a
%   3D unit cell / RVE with axis-aligned, rectangular faces. The routine:
%     1) Detects boundary nodes and classifies them into *opposite planes*,
%        their *intersection lines* (edges), and *corner points*.
%     2) Builds consistent **master/slave node pairs** across opposite faces.
%     3) Assembles the **constraint matrix G** such that, in incremental form,
%           d(u_slave) = G * d(u_master)
%        (rigid translation part dR is not used here; macro kinematics can be
%         applied externally if needed).
%     4) Reorders the global DOFs and returns a sparse matrix **A** that
%        enforces the constraints by substitution:
%           [u_slave; u_master; u_free]  ->  A * [u_master; u_free]
%
%   The function also defines a **reference configuration** (corner or center)
%   and returns `COORref`, i.e., coordinates measured w.r.t. that point, which
%   is useful to apply macroscopic strain-driven displacements externally as
%   u_macro(X) = E_macro * COORref.
%
% SCOPE & ASSUMPTIONS
%   - Geometry: Orthogonal, axis-aligned hexahedral RVE (faces parallel to
%     coordinate planes X=const, Y=const, Z=const). Non-skewed bounding box.
%   - Periodicity: Node-to-node pairing is assumed feasible (same surface
%     mesh on opposite faces up to translation).
%   - Large-strain framework: Intended for large-deformation solvers; this
%     routine only builds kinematic couplings (it is agnostic to the
%     constitutive law).
%
% INPUTS
%   COOR : (nnode × 3) nodal coordinates.
%   CN   : (nel × nnel) volume connectivity (used to infer a robust tolerance).
%   CNb  : (nel_b × nnel_b) boundary connectivity (to identify boundary nodes).
%   DATA : struct with optional fields:
%          - .TOL_deter_BOUNDARYNODES : numeric tolerance for classifying nodes
%            onto planes/lines/points. If empty/absent, `ChooseTolerance(CN,COOR)`
%            is used to infer a sensible value from mesh scales.
%          - .REFERENCE_POINT : 'CENTER' (default) or 'CORNER'.
%              'CENTER' -> reference at the geometric center of the RVE box.
%              'CORNER' -> reference at the first detected corner node.
%
% OUTPUTS
%   DISP_CONDITIONS : struct with fields
%       .A    : (ndof × ndofL) sparse substitution matrix enforcing PBCs by
%               mapping the *kept* DOFs [DOFm; DOFf] to the full set:
%                   u_full = A * u_kept,  where u_kept = [u_master; u_free].
%       .G    : (3*ns × 3*nm) block matrix such that d(u_slave)=G*d(u_master),
%               assembled by repeating a scalar coupling Gi on the 3 Cartesian
%               components.
%       .DOFr : vector of slave DOF indices (eliminated by substitution).
%       .DOFm : vector of master DOF indices (retained).
%       .DOFf : vector of free (interior) DOF indices (retained).
%       .DOFl : concatenation [DOFm; DOFf] defining the kept DOF order used by A.
%
%   COORref : (3 × nnode) coordinates shifted by the chosen reference point C:
%             COORref = (COOR' - C), with C = center or corner. Useful for
%             external application of macro kinematics.
%
%   DATA : input struct, augmented with:
%          - .VOL_RVE : RVE box volume = (xmax-xmin)*(ymax-ymin)*(zmax-zmin).
%          - .CENTER  : [xC yC zC], the geometric center of the box.
%
% METHOD (high-level steps)
%   1) Boundary box & reference:
%      - Extract boundary nodes from CNb and compute xmin/xmax, ymin/ymax, zmin/zmax.
%      - Store the box center and volume in DATA; choose reference point C.
%   2) Classification of boundary entities:
%      - Using tolerance TOL, classify boundary nodes into six planes
%        (X=min/max, Y=min/max, Z=min/max).
%      - Derive intersection *lines* (12 edges) and *corner* points (8).
%      - Remove overlaps so that sets are disjoint: points ⊂ lines ⊂ planes.
%   3) Master/Slave pairing:
%      - Call `MasterSlavesSets(...)` to create consistent pairings between
%        opposite planes and to assign which side is master vs. slave.
%   4) Coupling matrix:
%      - `GmatrixSlaveMaster(...)` builds the scalar coupling Gi from
%        (slave, master, COORref) pairings, then G is expanded to 3D by
%        block-diagonal replication on the components.
%   5) DOF sets and substitution matrix:
%      - Build DOFr (slaves), DOFm (masters), DOFf (free). Define kept order
%        DOFl = [DOFm; DOFf].
%      - Assemble sparse A so that:
%          A(DOFr, DOFm) = G,  A(DOFm, DOFm) = I,  A(DOFf, DOFf) = I,
%        and all other blocks are zero; then restrict columns to DOFl.
%
% HOW TO USE (typical pattern)
%   [BC, COORref, DATA] = Periodic_BoundaryCOND_LARGE(COOR, CN, CNb, DATA);
%   % Keep unknowns only for master+free DOFs:
%   Kred = BC.A' * K * BC.A;           % condensed stiffness
%   Fred = BC.A' * F;                  % condensed force
%   ured = solve(Kred, Fred);          % your solver
%   ufull = BC.A * ured;               % recover full displacement vector
%
%   % To drive with a macroscopic strain E_macro, build u_macro externally, e.g.:
%   % u_macro = kron(eye(3), E_macro) * COORref(:);  % or your preferred scheme
%   % and superpose or enforce as needed in your solver pipeline.
%
% DEPENDENCIES (must be available on MATLAB path)
%   DefaultField, ChooseTolerance,
%   DetermineePlanesPeriodic, DetermineLinesPeriodic_OLD, DeterminePointsPeriodic,
%   RemoveLinesPlanes, RemovePointsLine,
%   MasterSlavesSets, GmatrixSlaveMaster
%
% NUMERICAL NOTES / PITFALLS
%   - Tolerance TOL is critical. Too small → misclassification (missed pairs);
%     too large → wrong pairings across planes. If your mesh is coarse/irregular,
%     override DATA.TOL_deter_BOUNDARYNODES.
%   - Mesh compatibility: Opposite faces should be topologically compatible for
%     node-to-node pairing. Otherwise G construction can fail or produce
%     ill-posed constraints.
%   - Geometry: The algorithm assumes faces are aligned with global axes.
%     For skewed or curved periodic boundaries, a different pairing strategy
%     (e.g., projection-based or search via periodic images) is required.
%   - Large rotations: PBC kinematics here are linear couplings of nodal
%     displacements; they remain valid as long as your global formulation
%     consistently handles large-deformation measures elsewhere.
%
% OUTPUT CONSISTENCY CHECKS (suggested)
%   - size(BC.G) = [3*#slaveNodes, 3*#masterNodes]
%   - BC.A is ndof × (length(BC.DOFm)+length(BC.DOFf))
%   - All slave DOFs are eliminated: setdiff(1:ndof, BC.DOFl) == BC.DOFr
%
% HISTORY
%   Original idea (small strains & 2D/3D variants): Periodic_BoundaryCOND.m
%   This version adapted for large-deformation workflows and explicit A/G
%   assembly. Initial author: Joaquín A. Hernández, 29-Oct-2015.
%   This header updated: 05-Oct-2025.
%==========================================================================

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
DATA = DefaultField(DATA,'REFERENCE_POINT','CENTER') ; 
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

 
%  
% uB = MACRODEF*DIFFcoor';
% 
% % dispMACRO
% %dbstop('99')
% dispMACRO = MACRODEF*COORref ;
% dispMACRO = dispMACRO(:); 



DOFr = zeros(3*length(nods),1) ;
%dR = zeros(3*length(nods),1) ;
DOFm = zeros(3*length(nodm),1) ;
for i=1:3
    % dR = [dR; uB(i,:)'];
    ROWS = i:3:3*length(nods) ;
    DOFr(ROWS) = [  3*(nods-1)+i];
 %   dR(ROWS) = uB(i,:)' ;
    ROWS = i:3:3*length(nodm) ;
    DOFm(ROWS) = [ 3*(nodm-1)+i];
    
end

ndof = prod(size(COOR)) ; 
DOFf = setdiff(1:ndof,[DOFr;DOFm])' ; 

% Matrix A 
DOFl = [DOFm;DOFf] ; 
ndofL = length(DOFl) ; 
A = sparse(ndof,ndof) ;  % First we make it ndof,ndof

ROWS = DOFr ; 
%COLS = 1:length(DOFm) ; 
COLS = DOFm ;
A(ROWS,COLS ) = G; 
% 
ROWS = DOFm ; 
COLS = DOFm ;
%
A(ROWS,COLS ) = speye(length(ROWS)); 
% 
ROWS = DOFf ; 
COLS = DOFf ;  
A(ROWS,COLS ) = speye(length(ROWS)); 
A = A(:,DOFl) ; 

DISP_CONDITIONS.A = A; 
DISP_CONDITIONS.G = G; 

DISP_CONDITIONS.DOFr = DOFr; 
DISP_CONDITIONS.DOFl = DOFl; 
DISP_CONDITIONS.DOFm = DOFm; 
DISP_CONDITIONS.DOFf = DOFf; 











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