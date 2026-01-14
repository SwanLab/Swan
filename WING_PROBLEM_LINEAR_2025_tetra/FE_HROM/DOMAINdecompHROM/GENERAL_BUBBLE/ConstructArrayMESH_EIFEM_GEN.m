function MESH = ConstructArrayMESH_EIFEM_GEN(SizeArray,ScaleFactor,DEGREE_POLYNOMIAL,COOR,CN_lines,LINES_BOUNDARY)
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
% 0. Add dependencies and load temporary input (if function called without arguments)
% --------------------------------------------------------------------------------------------------

% --------------------------------------------------------------------------------------------------
% 1. Define corner geometry of a single auxetic unit cell (unscaled)
% --------------------------------------------------------------------------------------------------

% --------------------------------------------------------------------------------------------------
% 2. Interpolate edges of the polygon using specified polynomial order
%    Outputs nodal coordinates and connectivity along edges
% --------------------------------------------------------------------------------------------------

% --------------------------------------------------------------------------------------------------
% 3. Scale the interpolated coordinates by a user-defined factor
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
%    - ElementsBND_2: Top side (uses interpolation nodes)
%    - ElementsBND_3: Right side
%    - ElementsBND_4: Bottom side (uses interpolation nodes)
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
    load('tmp1.mat')
    SizeArray = [2,2] ; 
end
 
[COORquad, LINES] = generateInterpEdges_GEN2D(COOR,CN_lines, DEGREE_POLYNOMIAL) ;


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


% 
% PLOT_FIGURES = 0;
% 
% if PLOT_FIGURES == 1
%     
%     figure(1)
%     hold on
%     xlabel('x')
%     ylabel('y')
%     axis equal
%     
%     for inode = 1:size(COORcell,1)
%         xPLOT = [COORcell(inode,1)] ;
%         yPLOT = [COORcell(inode,2)] ;
%         plot(xPLOT,yPLOT)
%         text(xPLOT,yPLOT,num2str(inode) )
%         
%     end
%     
% end



[MESH.COOR,MESH.CN,ARRAY_CELLS] = MeshCELLrepeat2D(COORcell,TRANSLATIONS,SizeArray) ;
MESH.TypeElement = 'GIVEN_BY_USER' ;
MESH.MaterialType = ones(size(MESH.CN,1),1) ;

% AUXILIAR MESH FOR PLOTTING INTERFACES
MESH.AUXILIAR.CN  = [] ;
MESH.AUXILIAR.TypeElement = 'Linear' ;

%      Boundary connectivities:
%    - ElementsBND{1}: Left side
%    - ElementsBND{2}: Top side  
%    - ElementsBND{3}: Right side
%    - ElementsBND{4}: Bottom side  
nFACES = 4; % Number of faces of the mesh determined by  ARRAY_CELLS
ElementsBND =cell(nFACES,1) ; % Cell array containing all the line elements (determined by two nodes) of each face
% -------------------------------------------------------------------------------------------------------------------
line_ELEMENTS = cell(nFACES,1) ; 
line_ELEMENTS{1} = ARRAY_CELLS(:,1) ;  % These are the index of the EIF element that has edges in the boundary ibound = 1,2,3,4
line_ELEMENTS{2} = ARRAY_CELLS(end,:) ; 
line_ELEMENTS{3} =  ARRAY_CELLS(:,end) ; 
line_ELEMENTS{4} =  ARRAY_CELLS(1,:) ; 
% 
LocalElementsBoundary = cell(nFACES,1) ;  
for ifaces = 1:nFACES
LocalElementsBoundary{ifaces} = LINES(LINES_BOUNDARY{ifaces})  ; 
end
 


 

for iface= 1:nFACES
   % Loop over number of the boundaries whose nodes we wish to identify
   for iLINESloc = 1:length(LocalElementsBoundary{iface})
       % Loop over number of "sub-segments" of the face
       SEGMENT_loc = LocalElementsBoundary{iface}{iLINESloc} ;
       for inodes = 1:(length(SEGMENT_loc)-1)
           % Loop over number of nodes of the segment (2 for linear elements, 3 for quadratic...)
           SEGMENT_loc_2nodes = SEGMENT_loc([inodes,inodes+1]) ; 
           ElementsBND{iface} = [ ElementsBND{iface};MESH.CN(line_ELEMENTS{iface}, SEGMENT_loc_2nodes) ] ; 
           MESH.AUXILIAR.CN = [MESH.AUXILIAR.CN ; MESH.CN(:,SEGMENT_loc_2nodes)] ;
       end
   end
    
end

MESH.CNb = [ElementsBND{1}; ElementsBND{2};ElementsBND{3}; ElementsBND{4}] ;
MESH.TypeElementB = 'Linear';

MESH.NODES_FACES = {unique(ElementsBND{1}(:)),unique(ElementsBND{2}(:)),unique(ElementsBND{3}(:)),unique(ElementsBND{4}(:))} ;
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

% 
% if PLOT_FIGURES == 1
%     
%     figure(2)
%     hold on
%     xlabel('x')
%     ylabel('y')
%     axis equal
%     
%     for icoor = 1:size(COOR,1)
%         plot(COOR(icoor,1),COOR(icoor,2)) ;
%         text(COOR(icoor,1),COOR(icoor,2),num2str(icoor))
%     end
%     
% end
