function MESH = ConstructArrayMESH_cellGEN(COOR,SizeArray,ScaleFactor,ORDER_pol)
% Tailoring a  coarse-scale mesh for an AUXETIC CELL
% 4-APRIL-2025, FRIDAY GREEN'S ARIBAU, BARCELONA/11-mAY-2025, sUNDAY,
% bALMES 185, bARCELONA
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/107_EIFEM_plast3D/08_GENERAL_MODES.mlx
% /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/109_EIFEM_largeROT/12_REPETITIVE_TRAIN.mlx
% --------------------------------------------------------------------------------------------------------------------
% Comments by ChatGPT-4
% --------------------------------------------------------------------------------------------------------
% Constructs a 2D mesh by replicating a basic auxetic cell pattern.
% Inputs:
%   SizeArray   - [nX, nY] array indicating number of repetitions in X and Y.
%   ScaleFactor - Scaling factor to adjust the size of the basic cell.
%   ORDER_pol   - Polynomial interpolation order for each edge.
%
% Outputs:
%   MESH        - Struct containing mesh coordinates, connectivities, boundaries, and normals.

% --------------------------------------------------------------------------------------------------
% 1. COOR:   corner geometry of a single auxetic unit cell (unscaled)
% --------------------------------------------------------------------------------------------------

% --------------------------------------------------------------------------------------------------
% 2. Interpolate edges of the polygon using specified polynomial order
%    Outputs nodal coordinates and connectivity along edges (ORDER_pol)
% --------------------------------------------------------------------------------------------------
% NOTE: One pair of edges (specifically, edges 1 and 3) is assumed to be linear,
% --------------------------------------------------------------------------------------------------
% 3. Scale the interpolated coordinates by a user-defined factor (ScaleFactor)
% --------------------------------------------------------------------------------------------------

% --------------------------------------------------------------------------------------------------
% 4. Compute the horizontal and vertical translation vectors
%    These define how to tile the unit cell to form a larger mesh
% --------------------------------------------------------------------------------------------------

% --------------------------------------------------------------------------------------------------
% 5. Optional: Plot the geometry of the single cell for debugging/visualization
% --------------------------------------------------------------------------------------------------

% --------------------------------------------------------------------------------------------------
% 6. Replicate the scaled unit cell across the domain using translation vectors
%    The resulting mesh is stored in MESH.COOR and MESH.CN
% --------------------------------------------------------------------------------------------------

% --------------------------------------------------------------------------------------------------
% 7. Initialize element type and material identifiers
% --------------------------------------------------------------------------------------------------

% --------------------------------------------------------------------------------------------------
% 8. Initialize auxiliary connectivity matrix (used for visualizing edges/faces)
% --------------------------------------------------------------------------------------------------

% --------------------------------------------------------------------------------------------------
% 9. Construct boundary connectivities:
%    - ElementsBND_1: Left side
%    - ElementsBND_2: Bottom side (uses interpolation nodes)
%    - ElementsBND_3: Right side
%    - ElementsBND_4: Top side (uses interpolation nodes)
% --------------------------------------------------------------------------------------------------

% --------------------------------------------------------------------------------------------------
% 10. Store global boundary connectivities and associated element type
% --------------------------------------------------------------------------------------------------

% --------------------------------------------------------------------------------------------------
% 11. Organize boundary nodes per face and store in MESH.NODES_FACES
% --------------------------------------------------------------------------------------------------

% --------------------------------------------------------------------------------------------------
% 12. Identify which elements correspond to each face and store connectivity indices
% --------------------------------------------------------------------------------------------------

% --------------------------------------------------------------------------------------------------
% 13. Compute outward unit normal and tangent vectors for each boundary face
% --------------------------------------------------------------------------------------------------

% --------------------------------------------------------------------------------------------------
% 14. Optional: Plot original polygon (pre-scaling and pre-replication) for verification
% --------------------------------------------------------------------------------------------------
 
if ~exist('Quadratic3N')
    addpath(genpath(getenv('MATLAB_CODES_FEHROM')))  ;
end
if nargin == 0
    load('tmp.mat')
end


% 
% COOR = [  0.3780415430267064 -0.8434421364985166
%     2.6029673590504463 -0.8434421364985166
%     2.981008902077152 0.0
%     2.981008902077152 0.3145400593471811
%     2.6029673590504463 1.1579821958456977
%     0.3780415430267064 1.1579821958456977
%     0.0 0.3145400593471811
%     0.0 0.0  ] ;



[COORquad, LINES] = generateInterpolatedEdges(COOR, ORDER_pol) ;


COORcell = COORquad*ScaleFactor ;

% TRANSLATION VECTOR

idim   = 1;
[maxLOC,iii_max] = max(COORcell(:,idim)) ;
[minLOC,iii_min] = min(COORcell(:,idim)) ;
TRANSLATIONS{idim} = [maxLOC-minLOC,0];

idim   = 2;
[maxLOC,iii_max] = max(COORcell(:,idim)) ;
[minLOC,iii_min] = min(COORcell(:,idim)) ;
TRANSLATIONS{idim} =  -[0,maxLOC-minLOC];



PLOT_FIGURES = 0;

if PLOT_FIGURES == 1
    
    figure(1)
    hold on
    xlabel('x')
    ylabel('y')
    axis equal
    
    for inode = 1:size(COORcell,1)
        xPLOT = [COORcell(inode,1)] ;
        yPLOT = [COORcell(inode,2)] ;
        plot(xPLOT,yPLOT)
        text(xPLOT,yPLOT,num2str(inode) )
        
    end
    
end



[MESH.COOR,MESH.CN,ARRAY_CELLS] = MeshCELLrepeat2D(COORcell,TRANSLATIONS,SizeArray) ;
MESH.TypeElement = 'GIVEN_BY_USER' ;
MESH.MaterialType = ones(size(MESH.CN,1),1) ;

% AUXILIAR MESH FOR PLOTTING INTERFACES
%
% MESH.AUXILIAR.CN = [MESH.CN(:,[1,2]) ; MESH.CN(:,2:3); MESH.CN(:,3:4); MESH.CN(:,4:5); ...
%     MESH.CN(:,6:7); MESH.CN(:,8:9);  MESH.CN(:,9:10); MESH.CN(:,10:11); MESH.CN(:,11:12);...
%        MESH.CN(:,13:14)] ;

MESH.AUXILIAR.CN  = [] ;
MESH.AUXILIAR.TypeElement = 'Linear' ;

% BOUNDARY ELEMENTS % ----------------------------------------
% LINE 1
NNODE_coarse= size(COORcell,1) ;
line_1_ELEMENTS = ARRAY_CELLS(:,1) ;
LocalElementsBoundary_INI = [NNODE_coarse-1,NNODE_coarse] ;
ElementsBND_1 = MESH.CN(line_1_ELEMENTS,LocalElementsBoundary_INI) ;
%Indexes_faces_bnd_element_1 = 1:size(ElementsBND_1,1) ;
% LINE 2
line_2_ELEMENTS = ARRAY_CELLS(end,:) ;

ElementsBND_2 = [];

for i = 1:(length(LINES{1})-1)
    LocalElementsBoundary = [i, i+1];
    MESH.AUXILIAR.CN = [MESH.AUXILIAR.CN ; MESH.CN(:,LocalElementsBoundary)] ;
    ElementsBND_2 = [ElementsBND_2; MESH.CN(line_2_ELEMENTS, LocalElementsBoundary)];
end


% LINE 3
inodeLOC = LocalElementsBoundary(end)+1 ;
line_3_ELEMENTS = ARRAY_CELLS(:,end) ;
LocalElementsBoundary = [inodeLOC,inodeLOC+1] ;
ElementsBND_3 = MESH.CN(line_3_ELEMENTS,LocalElementsBoundary) ;
MESH.AUXILIAR.CN = [MESH.AUXILIAR.CN ; MESH.CN(:,LocalElementsBoundary)] ;

% LINE 4
inodeLOC = LocalElementsBoundary(end)+1 ;
line_4_ELEMENTS = ARRAY_CELLS(1,:) ;
ElementsBND_4 = [] ;
for i = inodeLOC:(inodeLOC + length(LINES{3})-2)
    LocalElementsBoundary = [i, i+1];
    ElementsBND_4 = [ElementsBND_4; MESH.CN(line_4_ELEMENTS, LocalElementsBoundary)];
    MESH.AUXILIAR.CN = [MESH.AUXILIAR.CN ; MESH.CN(:,LocalElementsBoundary)] ;
end

%
MESH.AUXILIAR.CN = [MESH.AUXILIAR.CN ; MESH.CN(:,LocalElementsBoundary_INI)] ;



MESH.CNb = [ElementsBND_1; ElementsBND_2;ElementsBND_3; ElementsBND_4] ;
MESH.TypeElementB = 'Linear';

MESH.NODES_FACES = {unique(ElementsBND_1(:)),unique(ElementsBND_2(:)),unique(ElementsBND_3(:)),unique(ElementsBND_4(:))} ;
MESH.NODES_LINES = MESH.NODES_FACES;


Indexes_faces_bnd_element = cell(1,length(MESH.NODES_FACES)) ;

for iface = 1:length(MESH.NODES_FACES)
    [dummy, setBelemLOC]= ElemBnd(MESH.CNb,MESH.NODES_FACES{iface}); % elements face "iface"
    % Connectivities faces f1 and f2
    CONNECTb_faces_iface = MESH.CNb(setBelemLOC,:) ;
    % The above is the connectivity matrix for the nodes of face "iface"
    Indexes_faces_bnd_element{iface} = setBelemLOC ;
    
    
    
    % KW:NORMALSlocal
    [NORMALv{iface},TANGENTv{iface}] = NormalsBoundaryLocal(MESH.COOR,CONNECTb_faces_iface )  ;
    
    
end
MESH.Indexes_faces_bnd_element = Indexes_faces_bnd_element;
MESH.NORMALv = NORMALv;
MESH.TANGENTv = TANGENTv;


if PLOT_FIGURES == 1
    
    figure(2)
    hold on
    xlabel('x')
    ylabel('y')
    axis equal
    
    for icoor = 1:size(COOR,1)
        plot(COOR(icoor,1),COOR(icoor,2)) ;
        text(COOR(icoor,1),COOR(icoor,2),num2str(icoor))
    end
    
end
