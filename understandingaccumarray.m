%% Understanding accumarray
% clc; clear;

% Tenim dofs, una matriu columna 768x1. El contingut de la matriu indica el
% DOF.
% L'assignacio a la matriu sparse es fa tal que (num_fila, dof) = valor per
% obtenir una matriu Bdof = 768 x 135 (ntot x ndofglob)

%% Sparse
% S = sparse(i,j,s,m,n,nzmax) uses vectors i, j, and s to generate an
%     m-by-n sparse matrix such that S(i(k),j(k)) = s(k)
% Generes la matriu sparse Bdof 768x135 i la sumes a la matriu original Bt
% també 768x135.

%% Accumarray
% B = accumarray(IND,DATA) sums groups of data by accumulating elements
%     of a vector DATA according to the groups specified in IND.

% Think of it this way: 
salesman   = [  1;   2;   4;   2;   4];
amountsold = [101; 102; 103; 104; 105];
A = accumarray(salesman, amountsold);

index = [1 2 1;
        2 1 2;
        2 3 2;
        2 1 2;
        2 3 2];
data = [101;
        102;
        103;
        104;
        105];
A = accumarray(index, data)