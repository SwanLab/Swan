%%%% TOPOLOGY OPTIMIZATION CODE BY ANDREAS RIETZ, FEBRUARY 1999 %%%
function TopOriginal(nelx,nely,volfrac,penal);
% INITIALIZE
x(1:nely,1:nelx) = volfrac; 
loop = 0; 
change = 1.;
% START ITERATION
for loop = 1:80 
loop = loop + 1;
xold = x;
% FE-ANALYSIS
[U]=FE(nelx,nely,x,penal); 
% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
[KE] = lk;
c = 0.;
for ely = 1:nely
for elx = 1:nelx
n1 = (nely+1)*(elx-1)+ely; 
n2 = (nely+1)* elx +ely;
Ue = U([2*n1-1;2*n1; 2*n2-1;2*n2; 2*n2+1;2*n2+2; 2*n1+1;2*n1+2],1);
c = c + x(ely,elx)^penal*Ue'*KE*Ue;
dc(ely,elx) = -penal*x(ely,elx)^(penal-1)*Ue'*KE*Ue;
end
end
% DESIGN UPDATE BY THE OPTIMALITY CRITERIA METHOD
[x] = MMA(nelx,nely,x,volfrac,dc); 
% PRINT RESULTS
change = max(max(abs(x-xold)));
disp([' It.: ' sprintf('%4i',loop) ' Obj.: ' sprintf('%10.4f',c) ...
' Vol.: ' sprintf('%6.3f',sum(sum(x))/(nelx*nely)) ...
' ch.: ' sprintf('%6.3f',change )])
% PLOT DENSITIES
figure(1); colormap(flipud(gray)); imagesc(x); axis equal; axis tight; 
drawnow
axis off; colorbar; figure(1); drawnow
end 
%%%%%%%%%% MMA UPDATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xnew]=MMA(nelx,nely,x,volfrac,dc)
xlow=0.001; xhigh=1;
L=x-0.1*(xhigh-xlow)*ones(nely,nelx);
high=(x-L).^2.*-dc./(xlow-L).^2;
low=(x-L).^2.*-dc./(xhigh-L).^2; 
l2 = min(min(high));
l1 = max(max(low));
for i=1:50
lmid = 0.5*(l2+l1);
xnew=max(xlow,min(xhigh,L+abs(x-L).*sqrt(-dc./lmid)));
if sum(sum(xnew)) - volfrac*nelx*nely > 0;
l2 = lmid;
else
l1 = lmid;
end
end
%%%%%%%%%% FE-ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [U]=FE(nelx,nely,x,penal)
[KE] = lk; 
K = sparse(2*(nelx+1)*(nely+1), 2*(nelx+1)*(nely+1));
F = sparse(2*(nely+1)*(nelx+1),1); U = sparse(2*(nely+1)*(nelx+1),1);
for ely = 1:nely
for elx = 1:nelx
n1 = (nely+1)*(elx-1)+ely; 
n2 = (nely+1)* elx +ely;
edof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2];
K(edof,edof) = K(edof,edof) + x(ely,elx)^penal*KE;
end
end
% DEFINE LOADS AND SUPPORTS
F(2*(nelx+1)*(nely+1),1) = -1;
fixeddofs = [1:2*(nely+1)];
alldofs = [1:2*(nely+1)*(nelx+1)];
freedofs = setdiff(alldofs,fixeddofs);
% SOLVING
U(freedofs,:) = K(freedofs,freedofs) \ F(freedofs,:); 
U(fixeddofs,:)= 0;

ux=reshape(U(1:2:2*(nely+1)*(nelx+1)-1), ...
(nely+1),(nelx+1));
uy=reshape(U(2:2:2*(nely+1)*(nelx+1)), ...
(nely+1),(nelx+1));

figure(2);
quiver(ux,uy);
title('Displacements');


%%%%%%%%%% ELEMENT STIFFNESS MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [KE]=lk
E = 1.; 
nu = 0.3;
k=[ 1/2-nu/6 1/8+nu/8 -1/4-nu/12 -1/8+3*nu/8 ... 
-1/4+nu/12 -1/8-nu/8 nu/6 1/8-3*nu/8];
KE = E/(1-nu^2)*[ k(1) k(2) k(3) k(4) k(5) k(6) k(7) k(8)
k(2) k(1) k(8) k(7) k(6) k(5) k(4) k(3)
k(3) k(8) k(1) k(6) k(7) k(4) k(5) k(2)
k(4) k(7) k(6) k(1) k(8) k(3) k(2) k(5)
k(5) k(6) k(7) k(8) k(1) k(2) k(3) k(4)
k(6) k(5) k(4) k(3) k(2) k(1) k(8) k(7)
k(7) k(4) k(5) k(2) k(3) k(8) k(1) k(6)
k(8) k(3) k(2) k(5) k(4) k(7) k(6) k(1)];
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code was written by Andreas Rietz, Department of Mechanics %
% at Linkoping University, Sweden %
% %
% The syntax of the code has been revised since february 1999 to be more %
% similar to the optimization code described in %
% "A 99 line topology optimization code written in Matlab" %
% by Ole Sigmund (To appear in Structural Optimization). %
% The code as well as a postscript version of the paper can be %
% downloaded from the web-site: http://www.topopt.dtu.dk %
% That code includes in addition a mesh independancy filter which makes %
% the result more reliable. %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%