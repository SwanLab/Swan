clc
clear all 
 load('DATAWS/Celasgen_lam10lay_128.msh.mat')

% ---------------------------------------------
% Compare with theoretical results (isotropic)
% --------------------------------------------
E = 230e3 ; 
nu = 0.25 ;
h = 0.16 ; 
A11 = E*h/(1-nu^2) ; 
A = [A11 nu*A11 0 ; nu*A11 A11 0 ; 0 0 (1-nu)/2*A11]; 

A11 = E*h^3/(1-nu^2)/12 ; 
D = [A11 nu*A11 0 ; nu*A11 A11 0 ; 0 0 (1-nu)/2*A11]; 

A11 = E*h/(1-nu^2) ; 
K = 1;
AS = K*(1-nu)/2*[A11 0 ; 0 A11];


 
 disp('---------------------------')
disp('Matrix A')
disp('----------------------------')
Acomp = Celas(1:3,1:3); 
disp(Acomp)
disp('---------------------------')
disp('Matrix B')
disp('----------------------------')
Bcomp = Celas(1:3,4:6); 
disp(Bcomp)
disp('---------------------------')
disp('Matrix D')
disp('----------------------------')
Dcomp = Celas(4:6,4:6); 
disp(Dcomp)
disp('---------------------------')
disp('Matrix AS')
disp('----------------------------')
AScomp = Celas(7:8,7:8); 
disp(AScomp)


%%%% D
disp('------------------------')%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Difference between A and Acomp (%)')

diffA = abs(A-Acomp)./abs(Acomp)*100 ;
disp(diffA)

disp('------------------------')%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Difference between D and Dcomp (%)')

diffD = abs(D-Dcomp)./abs(Dcomp)*100 ;
disp(diffD)

disp('------------------------')%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Difference between AS and AScomp (%)')

diffAS = abs(AS-AScomp)./abs(AScomp)*100 ;
disp(diffAS)

