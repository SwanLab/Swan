function [Yrot,SingYrot] = AlignmentSubspaces(Y,SingY,Z)
% Alignment of Subspaces 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y = bsxfun(@times,Y',SingY)' ;   %  Matrix we wish to align (multiplied by weighting factors
[Z,~,~ ]= SVDT(Z,0) ; % Reference matrix
A = Y'*Z ; 
[u,~,v] = SVDT(A,0) ; 
h = u*v' ; 
Yrot = Y*h ; 
[Yrot,SingYrot,~] = SVDT(Yrot,0) ; 

  