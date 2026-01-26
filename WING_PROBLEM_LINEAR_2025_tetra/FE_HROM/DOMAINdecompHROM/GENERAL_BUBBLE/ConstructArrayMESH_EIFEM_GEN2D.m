function MESH = ConstructArrayMESH_EIFEM_GEN2D(SizeArray,ScaleFactor,DEGREE_POLYNOMIAL,COOR,CN_lines,LINES_BOUNDARY)
% -------------------------------------------------------------------------
% ConstructArrayMESH_EIFEM_GEN2D
% -------------------------------------------------------------------------
% Constructs a repeated 2D mesh structure suitable for EIFEM (Empirical 
% Interscale Finite Element Method) simulations by tiling a reference 
% mesh cell across a structured array. Also identifies boundary edges, 
% nodes, and normal/tangent directions at the outer faces of the domain.
%
% INPUT:
%   - SizeArray:        [nRows, nCols] defining the number of repetitions
%                       of the unit cell in each spatial direction.
%   - ScaleFactor:      Scalar to scale the size of the reference geometry.
%   - DEGREE_POLYNOMIAL: Polynomial degree of interpolation for edges 
%                       (e.g., 1 = linear, 2 = quadratic).
%   - COOR:             Node coordinates of the reference mesh (Nx2).
%   - CN_lines:         Connectivity array (Mx2) defining segments/edges
%                       of the reference mesh.
%   - LINES_BOUNDARY:   Cell array of 4 entries with indices into CN_lines,
%                       identifying edges along:
%                         {1} Left, {2} Top, {3} Right, {4} Bottom.
%
% OUTPUT:
%   - MESH:             Structure containing the generated mesh fields:
%       - MESH.COOR:        Node coordinates of the full repeated mesh.
%       - MESH.CN:          Element connectivity matrix for each repeated cell.
%       - MESH.CNb:         Connectivity matrix of boundary edges.
%       - MESH.TypeElement: String describing the type of elements used.
%       - MESH.MaterialType: Material ID per element (default: 1).
%       - MESH.NODES_FACES: Cell array of node IDs for each of the 4 sides.
%       - MESH.NODES_LINES: Alias of NODES_FACES (may evolve).
%       - MESH.AUXILIAR:    Struct for auxiliary plotting/connectivity info.
%       - MESH.NORMALv:     Outward normal vectors (per boundary face).
%       - MESH.TANGENTv:    Tangent vectors (per boundary face).
%       - MESH.Indexes_faces_bnd_element: Indices of elements touching
%                                         each face, useful for interface methods.
%
% USAGE NOTES:
%   - Assumes periodic tiling of a base mesh cell over a 2D array.
%   - Normal and tangent vectors are useful for interface formulations.
%   - Used within the EIFEM framework to define structured microstructures.
%
% AUTHOR:
%   Joaquín A. Hernández, UPC, Campus Terrassa
%   27-May-2025
%   Commented by ChatGPT4
% -------------------------------------------------------------------------

% --------------------------------------------------------------------------------------------------
 if ~exist('Quadratic3N')
    addpath(genpath(getenv('MATLAB_CODES_FEHROM')))  ;
end
if nargin == 0
    load('tmp1.mat')
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
