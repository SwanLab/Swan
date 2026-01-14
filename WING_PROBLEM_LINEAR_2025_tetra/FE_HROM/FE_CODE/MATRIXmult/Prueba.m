clc
clear all 

% Example of how to multiply a matrix using bsxfun, permute and sum 
% ----------------------------------------------------------------------
A = [5 6; 5 8] ;  % m x n
B = ones(2,5) ;   % n x r 
C = A*B   ;  %   m x r 
% and the 1st dimension of the second matrix. 
Aperm = permute(A,[1,3,2])  ;  % We move the second dimension to the 3rd dimension  --> m x 1 x n  
Bperm = permute(B,[3,2,1])   ;   %  First dimension to the 3rd dimension          --->  1 x r x n
Cperm = bsxfun(@times,Aperm,Bperm)                                                % ---> m x r x n
Cnew = sum(Cperm,3)  ;                                                            % ---> m x r 