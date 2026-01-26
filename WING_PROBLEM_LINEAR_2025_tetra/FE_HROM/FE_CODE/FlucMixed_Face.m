function [DOF,r,l,J] =   FlucMixed_Face(nodesf,U,...
    R,M)
ndim = 3;
DOF = small2large(nodesf,ndim) ;   % DOFS face 
% Product by mass matrix 
% ----------------------
% Ubar = zeros(size(U)) ;
% Uout = zeros(size(Uin)) ;
% Rbar = zeros(size(R)) ;
% for idim =1:ndim
%     INDLOC =idim:ndim:size(U,1) ;
%     Ubar(INDLOC,:) = M*U(INDLOC,:) ;
%     Uout(INDLOC,:) = M*Uin(INDLOC,:) ;
%     Rbar(INDLOC,:) = M*R(INDLOC,:) ;
% end
Ubar = M*U  ; 

% Matrix Q
coeff = (Ubar'*U)\Ubar';
Q = U*coeff ;
Q = speye(size(Q))-Q ;
%Select p linearly independent rows from Ubar

[~,l]=licols([Ubar ]') ; %
r = setdiff(1:length(DOF),l) ; 
%
J =  Q(r,r)\Q(r,l) ;
