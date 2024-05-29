% Checking double contraction
clc; clear

A = reshape(9:-1:1, 3,3);
B = reshape(1:81, 3,3,3,3);
C = reshape(1:9, 3,3);

parc1 = einsum(A,B, 'ab,abkl->kl');
% parc2 = einsum(parc1, C,'kl,kl')
parc2 = einsum(A,B,C,'ab,abkl,kl')