function [Nshape,COORnodes,OTHER_OUTPUT] =  LinearPieceWise2Dmodes_EIFEM(INFO_INTERFACE_MESH,CENTROID,COORbnd)
% JAHO, 2-APRIL-2025, UPC, CAMPUS NORD, BARCELONA
% Computation of Nshape for the interface ficti. modes, in EIFEM
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/107_EIFEM_plast3D/08_GENERAL_MODES.mlx
% --------------------------------------------------------------------------------------
if nargin == 0
    load('tmp.mat')
end
OTHER_OUTPUT = [] ; 

COORnodes = INFO_INTERFACE_MESH.COOR; %(:,2:end) ;
ElemBnd = INFO_INTERFACE_MESH.LINES ;

nmodes = size(COORnodes,1) ; % Number of modes
npoints = size(COORbnd,1) ;

Nshape = zeros(npoints,nmodes) ;

% cOORDINATES REFERRED TO THE CENTROID
ndim = size(CENTROID,2);
for idim = 1:ndim
    COORnodes(:,idim) = COORnodes(:,idim) - CENTROID(idim) ;
end


for ielem = 1:size(ElemBnd,1)
    % We have to compute 2 shape functions per element
    % x = x_0 + lambda*uDIR
    iii = ElemBnd(ielem,1) ;
%    aa =  find(iii == LIST_NODES ) ;
    x1 = COORnodes(iii,:) ;
    iii = ElemBnd(ielem,2) ;
 %   aa =  find(iii == LIST_NODES ) ;
    x2 = COORnodes(iii,:) ;
    norm_x1x2 = norm(x2-x1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    
    % Compute the coordinates of all the boundary nodes with
    % respect to the first point of the interface we are considering
    COORbnd_rel = bsxfun(@minus,COORbnd',x1(:))' ;
    % Orthogonalize the vectors
    norm_COORbnd_rel = sqrt(sum(COORbnd_rel.^2,2)) ;
    % See if any of the points is the reference point
    TOL_zero  = 1e-8;
    [DDD,IIII]=min(norm_COORbnd_rel) ;
    CandidatePoints = 1:length(norm_COORbnd_rel) ;
    if DDD/norm_x1x2  <TOL_zero
        IndPointsInclude= IIII ;
        CandidatePoints = setdiff(1:length(norm_COORbnd_rel),IndPointsInclude) ;
    end
    
    
    
    uBND = bsxfun(@times,COORbnd_rel(CandidatePoints,:),1./norm_COORbnd_rel(CandidatePoints)) ;
    % Compute scalar product with  the vector that goes from x1 to x2
    u  = (x2-x1)/(norm(x2-x1)) ;
    PROY = u*uBND' ;
    % A necessary condition ---but not sufficient--- is that the nodes should
    % lie on the  line defined by x2 and x1
    TOL_zero = 1e-4;
    IndElements_1 = find(abs(PROY-1)/norm_x1x2<=TOL_zero );
    % The other condition is that COORbnd_rel divided by the norm of x2-x1
    % should be less than one
    Proy_u = u*COORbnd_rel(CandidatePoints(IndElements_1),:)'/norm_x1x2 ;
    TOL  = 1e-5 ;
    SubIndElements_2 = find(Proy_u <= (1+TOL) & Proy_u>=TOL); 
    IndPoints = IndElements_1(SubIndElements_2) ;
    
    IndPoints= unique([IndPointsInclude,CandidatePoints(IndPoints)]) ;
    
    % The   shape function 2*(ielem-1)+1 
    imode= 2*(ielem-1)+1; 
    % is obtained as (xP2-xP)/(x2P-x1P)
    % where xP is the projection onto u of the boundary nodes of the
    % segment 
    Nimode = -bsxfun(@minus,COORbnd(IndPoints,:)',x2(:))' ; 
    Nimode = u*Nimode'/norm_x1x2 ; 
   Nshape(IndPoints,imode)  = Nimode ; 
    
    % The   shape function 2*ielem 
    imode= 2*ielem ;  
    % is obtained as (x-x1P)/(x2P-x1P)
    % where xP is the projection onto u of the boundary nodes of the
    % segment 
    Nimode = bsxfun(@minus,COORbnd(IndPoints,:)',x1(:))' ; 
    Nimode = u*Nimode'/norm_x1x2 ; 
    Nshape(IndPoints,imode)  = Nimode ; 
    
    
end

