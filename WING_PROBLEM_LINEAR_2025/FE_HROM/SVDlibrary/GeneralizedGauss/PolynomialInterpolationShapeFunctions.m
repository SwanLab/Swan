function [N,B] = PolynomialInterpolationShapeFunctions(COORelement,n,COORevaluate)
% Given COORelement -> nnode x ndim, this  function returns the  coefficients of the
% polynomial shape function of order n(1)x n(2), such that
% N_inode(COORelement(inode,1),COORelement(inode,2)) = 1 .
% Then, it evaluates the  shape functions at point COORevaluate --> N(COORevaluate)
% -------------------------------------------------------------

% p(x,y) = [x^0*y^0,x^1*y^0,x^2*y^0,  ...
%     x^0*y^1,x^1*y^1,x^2*y^1  ...
%    x^0*y^2,x^1*y^2,x^2*y^2    ]_{1x9} * A_{9,1}  =
%    P_i(x,y)*A
% P(x,y)=   [1,x,x^2,  ...
%     y,x*y,x^2*y  ...
%    y^2,x*y^2,x^2*y^2    ]_{1x9}
% JAHO, April-2020
% ----------------------------------------------------------------------------------
if nargin == 0
    load('tmp1.mat')
    COORevaluate = COORelement(1:2,:) ; 
end

ndim = size(COORelement,2) ;
nnode = size(COORelement,1) ;
P = CoordinateMatrixPolynomial(COORelement,n)  ;
COEFFSpol = P\eye(nnode);  % Coefficients of the polynomials
% Evaluating the shape functions a point COORevaluate
[Pevaluate PevaluateDER]= CoordinateMatrixPolynomial(COORevaluate,n)  ;
N = Pevaluate*COEFFSpol ;  % Shape functions at the given points COORevaluate 
B = cell(length(PevaluateDER),1) ; 
for idim = 1:length(PevaluateDER)
   B{idim} =   PevaluateDER{idim}*COEFFSpol ; 
end
 