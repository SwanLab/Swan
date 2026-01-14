function  [MESH] = ConstructArrayMESH_2D(SizeArray,MESHrve) ;
 
% -------------------------------------------------------------------------

% --------------------------------------------------------------------------------------------------
if ~exist('Quadratic3N')
    addpath(genpath(getenv('MATLAB_CODES_FEHROM')))  ;
end
if nargin == 0
    load('tmp.mat')
    SizeArray = [1,2] ;
end


%
px = SizeArray(2);
py = SizeArray(1);
MAT = MESHrve.MaterialType;

% Mesh repetition (by ChatGPT)
% [COOR_ALL, CN_ALL, MAT_ALL, DOMAINS_INDEXES] = ...
%     createPeriodicMesh2D_vectorized(MESHrve.COOR, MESHrve.CN, px, py, MAT);

VECTORIZED_VERSION =1 ; 

if VECTORIZED_VERSION == 1
    [COOR_ALL, CN_ALL,CNb_ALL MAT_ALL, DOMAINS_INDEXES] = ...
    createPeriodicMesh2D_withBoundary_CGPT_vect(MESHrve.COOR, MESHrve.CN, MESHrve.CNb, px, py, MAT);
else

 [COOR_ALL, CN_ALL,CNb_ALL MAT_ALL, DOMAINS_INDEXES] = ...
    createPeriodicMesh2D_withBoundary(MESHrve.COOR, MESHrve.CN, MESHrve.CNb, px, py, MAT);

end

MESH = MESHrve;
MESH.COOR = COOR_ALL ;
MESH.CN   = CN_ALL ;
MESH.CNb   = CNb_ALL ;

MESH.MaterialType = MAT_ALL ;




% BOUNDARIES
% BY DEFAULT, WE ARE GOING TO TAKE 4 BOUNDARIES FOR IMPOSING BOUNDARY
% CONDITIONS
% RIGHT (1), BOTTOM (2), LEFT(3), TOP (4)


BoundaryNodes= unique(CNb_ALL) ; 
% Extract X and Y
X = COOR_ALL(BoundaryNodes,1);
Y = COOR_ALL(BoundaryNodes,2);

% Extremes
x_min = min(X);
x_max = max(X);
y_min = min(Y);
y_max = max(Y);

tol = 1e-4*(x_max-x_min);


% Find nodes on each boundary line
nodes_left   = find(abs(X - x_min) < tol);
nodes_right  = find(abs(X - x_max) < tol);
nodes_bottom = find(abs(Y - y_min) < tol);
nodes_top    = find(abs(Y - y_max) < tol);

NODES_LINES = cell(1,8);
NODES_LINES{1} = BoundaryNodes(nodes_left) ;
NODES_LINES{2} = BoundaryNodes(nodes_bottom) ;
NODES_LINES{3} = BoundaryNodes(nodes_right) ;
NODES_LINES{4} = BoundaryNodes(nodes_top) ;

% These are the boundaries of the reference domain
for i=1:4
    iLOC = 4+i ;
    NODES_LINES{iLOC} = MESHrve.NODES_LINES{i};
end

% NODES_FACES
NODES_FACES = cell(1,length(DOMAINS_INDEXES)) ;

for idomains = 1:length(DOMAINS_INDEXES)
    ielemLOC = DOMAINS_INDEXES{idomains} ;
    NODES_FACES{idomains}  = unique(MESH.CN(ielemLOC,:)) ;
    
end

MESH.NODES_LINES = NODES_LINES;
MESH.NODES_FACES = NODES_FACES;


