function [DOFr,dR,OTHEROUTPUT] = MIXEDCONDtime_PLATESfromFEshape(DIRICHLET,DATA,ndim,MESH,GEOproperties,DATALOC)
% Goal. Determine DOFr and    dR(t) , as well as traction boundary
% conditions
% Case boundary conditions expressed as constrained boundary shape
% functions (PLATEs)
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/103_EIFEM_PLATES/02_CONV_SIZE_TRAC.mlx
%
% JAHO- 23-AUG-2O23
%--------------------------------------------------------------------------------------
if nargin == 0
    load('tmp.mat')
elseif nargin == 5
    DATALOC = [] ;
end
OTHEROUTPUT = [] ;

FACES_BND  = DATALOC.LabelEntitiesDefiningBoundary ; % FACES DEFINING THE BOUNDARIES OF THE FINE MESH DOMAIN
%
% These are the "corner" lines, the intersection between the lateral
% surfaces of the plate element
LINES_INDEX_FACE = {[FACES_BND(4),FACES_BND(1)],[FACES_BND(1),FACES_BND(2)],[FACES_BND(2),FACES_BND(3)],[FACES_BND(3),FACES_BND(4)]} ;


NODES_LINES =  cell(size(LINES_INDEX_FACE)) ;
NODES_CORNERS = cell(size(NODES_LINES)) ;
DeltaZ = zeros(length(NODES_LINES),1) ;

for ilines = 1:length(LINES_INDEX_FACE)
    iface = LINES_INDEX_FACE{ilines}(1) ;
    jface = LINES_INDEX_FACE{ilines}(2) ;
    % THESE ARE THE NODES OF THE EDGES OF THE TRAINING DOMAIN  (LINES)
    NODES_LINES{ilines} =    intersect(MESH.NODES_FACES{iface},MESH.NODES_FACES{jface}) ;
    
    % FROM THESE LINES, WE TAKE THE TOP AND BOTTOM NODES (CORNERS). IT IS
    % TACITLY ASSUMED THAT THE TRAINING DOMAIN IS A RECTANGULAR CUBOID,
    % WITH ITS THICKNESS ALONG Z AXIS SIGNIFICANTLY SMALLER THAN IN THE
    % OTHER TWO DIRECTIONS
    [zmax,cornerlocMAX] =  max(MESH.COOR(NODES_LINES{ilines},3)) ;
    [zmin,cornerlocMIN] =  min(MESH.COOR(NODES_LINES{ilines},3)) ;
    DeltaZ(ilines)= zmax-zmin ;
    NODES_CORNERS{ilines} = [NODES_LINES{ilines}(cornerlocMIN),NODES_LINES{ilines}(cornerlocMAX)]' ;
    
end


NODES_CORNERS = cell2mat(NODES_CORNERS) ;
NODES_CORNERS = NODES_CORNERS' ;
NODES_CORNERS = NODES_CORNERS(:) ;

NODES_FACES = MESH.NODES_FACES(FACES_BND) ;  %



NODES_FACES = NODES_FACES(:) ;               %
rnodLOC = unique(cell2mat(NODES_FACES)) ;  % LIST OF NODES WITH PRESCRIBED DISPLACEMENTS (BOUNDARY)
% ----------------------------------------------------------
COORbnd = MESH.COOR(rnodLOC,:) ;  % COORDINATES

% SHAPE FUNCTIONS
% ---------------
% See example in /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/ContinuousEmpiricalCubatureM/Paper_hernandez2021ecm/09_ASSESSMENTpaper/Fun3D/Dpoly3D_P5.m
%
MESHcoar_COOR = MESH.COOR(NODES_CORNERS,:) ;

[nnodeBND, ndim] = size(MESHcoar_COOR) ;
% Order of polynomial
nnodes_DIM = (nnodeBND)^(1/ndim)  ;  % Number of nodes per dimension
ORDER_POLYNOMIALS = (nnodes_DIM-1)*ones(1,ndim) ;
DATAshape = ShapeFunCoefficients(MESHcoar_COOR,ORDER_POLYNOMIALS) ;
DATAlocSHAPE.DATAshape  = DATAshape;
xLIM = [] ;
DATAlocSHAPE.ORDER_POLYNOMIALS = ORDER_POLYNOMIALS;
[Nshape, ~,~ ]=    ShapeFunctionFE(xLIM,COORbnd,DATAlocSHAPE) ;

% NEXT STEP: DETERMINE THE "MODES" ASSOCIATED TO THESE SHAPE FUNCTIONS (8x 3 = 24 )
Vfe= zeros(size(Nshape,1)*ndim,size(Nshape,2)*ndim) ;
for idim = 1:ndim
    Vfe(idim:ndim:end,idim:ndim:end) = Nshape  ;
end
% We need to construct a  24 x 20 matrix relating Vfe with the desired
% modes U
% U = Vfe*A
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/101_MULTI2D_2023/02_PLATES/README_PLATES.mlx
Vplate = TransformHEXA_COORD_2_PLATE_COOR(Vfe,DeltaZ) ;
U = Vplate*DIRICHLET.AMPLITUDE ;
% Displacement pattern
DOFr = small2large(rnodLOC,ndim) ;  % Prescribed DOFS
% U = zeros(length(DOFr),1) ;   % Space pattern
% for idim = 1:ndim
%    U(idim:ndim:end) = Nshape*DIRICHLET.AMPLITUDE(idim:ndim:end) ;
% end
a = DIRICHLET.TIMEFUN(DATA.STEPS) ;  % Time pattern
dR.U = U ;
dR.a = a ;

% The foregoing are but the same steps followed in
% DirichletCONDtime_PLATESfromFEshape.m to derive the
% MULTIPLICATION OF U BY THE MASS MATRIX TO DETERMINE  equivalent nodal
% forces
% We may follow the steps outlined in
% /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/LARGE_STRAINS/InputDataFunctions/PointForceGeneralized.m

% Loop over FACES_BND
NNODES = size(MESH.COOR,1) ;
NODAL_FORCES_TRACTIONS = zeros(ndim*NNODES,1) ;  % This is the vector of nodal forces due to the surface traction
PRESCRIBED_DISP = {} ; 
PRESCRIBED_DOFS = {} ; 

for ifaceLOC = 1:length(FACES_BND) % Loop over the 4 faces of the plate element
    iface = FACES_BND(ifaceLOC) ;
    %
    nodesf  = MESH.NODES_FACES{iface} ;  % These are the nodes of this face, global numbering 
    
    [~,INDLOCNODESbnd] = intersect(rnodLOC,nodesf) ; %  Nodes of the face (local, with respet to the list of boundary nodes)
    DOFf_local = small2large(INDLOCNODESbnd,ndim) ;   % These are the corresponding DOFs (local)
    LOCAL_DISPLACEMENT = U(DOFf_local) ;  % These are the nodal displacements
    % MAss matrix
    Mst = MESH.PROPERTIES_FACES{iface}.GeometricMassMatrix  ;
    M = sparse(size(Mst,1)*ndim,size(Mst,2)*ndim) ;
    for idim = 1:ndim
        M(idim:ndim:end,idim:ndim:end) = Mst ;
    end
    % NODAL FORCES
    LOCAL_FORCES= M*LOCAL_DISPLACEMENT ;
    
    %% ASSEMBLY
    DOFf = small2large(nodesf,ndim) ;   % These are the corresponding DOFs, in global numbering
    
    NODAL_FORCES_TRACTIONS(DOFf) =  NODAL_FORCES_TRACTIONS(DOFf) +LOCAL_FORCES ;
    
    % DIRICHLET BOUNDARY CONDITIONS 
    % 
    if norm(LOCAL_DISPLACEMENT) < 1e-10 
        % This is a fixed face. Accordingly, we set 
        PRESCRIBED_DOFS{end+1} = DOFf; 
        PRESCRIBED_DISP{end+1} = LOCAL_DISPLACEMENT ; 
    end
    
    
    
    
    
    
end

OTHEROUTPUT.Ftrac.U = NODAL_FORCES_TRACTIONS; 
OTHEROUTPUT.Ftrac.a = dR.a; 

% PRESCRIBED DISPLACEMENTS
% ------------------------
PRESCRIBED_DOFS = cell2mat(PRESCRIBED_DOFS')  ; 
PRESCRIBED_DISP = cell2mat(PRESCRIBED_DISP')  ; 
[DOFr,III] = unique(PRESCRIBED_DOFS) ; 
dR.U = PRESCRIBED_DISP(III) ; 


% DOFr = small2large(rnodLOC,ndim) ;  % Prescribed DOFS
% % U = zeros(length(DOFr),1) ;   % Space pattern
% % for idim = 1:ndim
% %    U(idim:ndim:end) = Nshape*DIRICHLET.AMPLITUDE(idim:ndim:end) ;
% % end
% a = DIRICHLET.TIMEFUN(DATA.STEPS) ;  % Time pattern
% dR.U = U ;
% dR.a = a ;



