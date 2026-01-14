function MESH = ConstructArrayMESH_4sides_2D(SizeArray,ScaleFactor,DATALOC)
% This function construct a FE mesh by repeating a single unit cell, in 2D
% 
%  
% 2-MAY-2025, FABRA BY SECRETS, FRIDAY. 
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/109_EIFEM_largeROT/12_REPETITIVE_TRAIN.mlx
% AND e.g. /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/109_EIFEM_largeROT/12_REPETITIVE_TRAIN/FE_nonSM2_rep.m
% --------------------------------------------------------------------------------------------------------------------
if ~exist('Quadratic3N')
    addpath(genpath(getenv('MATLAB_CODES_FEHROM')))  ; 
end
if nargin == 0 
   load('tmp.mat')
   NUMBER_OF_CELLS = [2,1];  
end

% LET US READ THE MESH GIVEN IN
% DATALOC.
DATALOC.RenumberElementsForEficiency = 0;  
MESH_1cell = GeometryMesh(DATALOC) ;
%    COOR: [1548×2 double]
%                              CN: [344×9 double]
%                     TypeElement: 'Quadrilateral'
%                             CNb: [250×3 double]
%                    TypeElementB: 'Linear'
%                    MaterialType: [344×1 double]
%                     NODES_FACES: {[9×1 double]  [45×1 double]  [9×1 double]  [45×1 double]}
%                     NODES_LINES: {[9×1 double]  [45×1 double]  [9×1 double]  [45×1 double]}
%                       NODES_DOM: {}
%       Indexes_faces_bnd_element: {[4×1 double]  [22×1 double]  [4×1 double]  [22×1 double]}
%      NormalBoundaryElementsFace: {[2×4 double]  [2×22 double]  [2×4 double]  [2×22 double]}
%       TangentBoundaryElementsFace: {[2×4 double]  [2×22 double]  [2×4 double]  [2×22 double]}
%       IndicesRenumberingElements: []
 
COORcell = MESH_1cell.COOR*ScaleFactor ; 
NodeLeftSide = MESH_1cell.NODES_LINES{1}(1);
NodeRighttSide = MESH_1cell.NODES_LINES{3}(1);
NodeTopSide =  MESH_1cell.NODES_LINES{4}(1);
NodeBottomSide =  MESH_1cell.NODES_LINES{2}(1);

% Translation in DIRECTION 1 
dx = COORcell(NodeRighttSide,1) - COORcell(NodeLeftSide,1)  ; 
TRANSLATIONS{1} = [dx,0];
dy = COORcell(NodeTopSide,1) - COORcell(NodeBottomSide,1)  ; 
TRANSLATIONS{2} = [0,-dy];

% Now we have to construct the final mesh by making tiled copies of the
% firs
 



[MESH.COOR,MESH.CN,ARRAY_CELLS] = MeshCELLrepeat2D(COORcell,TRANSLATIONS,SizeArray) ;
MESH.TypeElement = 'GIVEN_BY_USER' ; 
MESH.MaterialType = ones(size(MESH.CN,1),1) ; 
 
% BOUNDARY ELEMENTS % ----------------------------------------
% LINE 1 
line_1_ELEMENTS = ARRAY_CELLS(:,1) ; 
LocalElementsBoundary = [9,10] ; 
ElementsBND_1 = MESH.CN(line_1_ELEMENTS,LocalElementsBoundary) ; 
%Indexes_faces_bnd_element_1 = 1:size(ElementsBND_1,1) ; 
% LINE 2 
line_2_ELEMENTS = ARRAY_CELLS(end,:) ; 
LocalElementsBoundary = [1,2] ; 
ElementsBND_2 = MESH.CN(line_2_ELEMENTS,LocalElementsBoundary) ; 
LocalElementsBoundary = [2,3] ; 
ElementsBND_2 = [ElementsBND_2;MESH.CN(line_2_ELEMENTS,LocalElementsBoundary)] ; 
%Indexes_faces_bnd_element_2 = (1:size(ElementsBND_2,1)) +Indexes_faces_bnd_element_1(end)  ; 

% LINE 3 
line_3_ELEMENTS = ARRAY_CELLS(:,end) ; 
LocalElementsBoundary = [4,5] ; 
ElementsBND_3 = MESH.CN(line_3_ELEMENTS,LocalElementsBoundary) ;
%Indexes_faces_bnd_element_3 =( 1:size(ElementsBND_3,1) )+Indexes_faces_bnd_element_2(end)  ; 

% LINE 4 
line_4_ELEMENTS = ARRAY_CELLS(1,:) ; 
LocalElementsBoundary = [6,7] ; 
ElementsBND_4 = MESH.CN(line_4_ELEMENTS,LocalElementsBoundary) ; 
LocalElementsBoundary = [7,8] ; 
ElementsBND_4 = [ElementsBND_4;MESH.CN(line_4_ELEMENTS,LocalElementsBoundary)] ; 
%Indexes_faces_bnd_element_4 = (1:size(ElementsBND_4,1)) +Indexes_faces_bnd_element_3(end)  ; 


MESH.CNb = [ElementsBND_1; ElementsBND_2;ElementsBND_3; ElementsBND_4] ; 
MESH.TypeElementB = 'Linear'; 

MESH.NODES_FACES = {unique(ElementsBND_1(:)),unique(ElementsBND_2(:)),unique(ElementsBND_3(:)),unique(ElementsBND_4(:))} ;
MESH.NODES_LINES = MESH.NODES_FACES; 

%MESH.Indexes_faces_bnd_element = {Indexes_faces_bnd_element_1,Indexes_faces_bnd_element_2,...
 %   Indexes_faces_bnd_element_3,Indexes_faces_bnd_element_4} ; 

% AUXILIAR MESH FOR PLOTTING INTERFACES 

MESH.AUXILIAR.CN = [MESH.CN(:,[1,2]) ; MESH.CN(:,2:3); MESH.CN(:,4:5); MESH.CN(:,6:7);  MESH.CN(:,7:8); MESH.CN(:,9:10);  ] ; 
MESH.AUXILIAR.TypeElement = 'Linear' ; 

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

% BOUNDARY ELEMENTS 
% IN THIS MESH, 