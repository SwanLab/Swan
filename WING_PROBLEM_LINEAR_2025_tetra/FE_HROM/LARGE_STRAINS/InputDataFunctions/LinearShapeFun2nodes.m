function  [Nshape,xi,IndPoints] =   LinearShapeFun2nodes(ElemBnd_loc,COORnodes,COORbnd,Nshape)


x1 = COORnodes( ElemBnd_loc(1),:) ;
x2 = COORnodes(ElemBnd_loc(2),:) ;
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
% Normalization
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

% Normalized coordinates
xREF = (x1+x2)/2; 
COORbnd_rel_0 = 2*bsxfun(@minus,COORbnd',xREF(:))'/norm_x1x2 ;
xi = u*COORbnd_rel_0(IndPoints,:)' ;



% The   shape function 2*(ielem-1)+1
imode= ElemBnd_loc(1);
% is obtained as (xP2-xP)/(x2P-x1P)
% where xP is the projection onto u of the boundary nodes of the
% segment
Nimode = -bsxfun(@minus,COORbnd(IndPoints,:)',x2(:))' ;
Nimode = u*Nimode'/norm_x1x2 ;
Nshape(IndPoints,imode)  = Nimode ;

% The   shape function 2*ielem
imode=ElemBnd_loc(2) ;
% is obtained as (x-x1P)/(x2P-x1P)
% where xP is the projection onto u of the boundary nodes of the
% segment
Nimode = bsxfun(@minus,COORbnd(IndPoints,:)',x1(:))' ;
Nimode = u*Nimode'/norm_x1x2 ;
Nshape(IndPoints,imode)  = Nimode ;

 

