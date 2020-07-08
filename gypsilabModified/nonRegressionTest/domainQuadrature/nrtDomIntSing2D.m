%+========================================================================+
%|                                                                        |
%|            This script uses the GYPSILAB toolbox for Matlab            |
%|                                                                        |
%| COPYRIGHT : Matthieu Aussal (c) 2017-2018.                             |
%| PROPERTY  : Centre de Mathematiques Appliquees, Ecole polytechnique,   |
%| route de Saclay, 91128 Palaiseau, France. All rights reserved.         |
%| LICENCE   : This program is free software, distributed in the hope that|
%| it will be useful, but WITHOUT ANY WARRANTY. Natively, you can use,    |
%| redistribute and/or modify it under the terms of the GNU General Public|
%| License, as published by the Free Software Foundation (version 3 or    |
%| later,  http://www.gnu.org/licenses). For private use, dual licencing  |
%| is available, please contact us to activate a "pay for remove" option. |
%| CONTACT   : matthieu.aussal@polytechnique.edu                          |
%|             francois.alouges@polytechnique.edu                         |
%| WEBSITE   : www.cmap.polytechnique.fr/~aussal/gypsilab                 |
%|                                                                        |
%| Please acknowledge the gypsilab toolbox in programs or publications in |
%| which you use it.                                                      |
%|________________________________________________________________________|
%|   '&`   |                                                              |
%|    #    |   FILE       : nrtDomIntSing2D.m                             |
%|    #    |   VERSION    : 0.50                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 25.11.2018                                    |
%|  / 0 \  |   LAST MODIF : 25.11.2018                                    |
%| ( === ) |   SYNOPSIS   : Singular integration over a triangle          |
%|  `---'  |                                                              |
%+========================================================================+

% Cleaning
clear all
close all
clc

% Gypsilab path
run('../../addpathGypsilab.m')

% Edge definition
a = 0;
b = 1;
N = 100;

% Reference mesh (edge)
vtx  = [a 0 0 ; b 0 0];
elt  = [1 2];
mesh = msh(vtx,elt);

% Observation points (100 random & 11 segments points)
X = [-2 + 5*rand(N,3);
    [(-2:0.2:3)' zeros(26,2)] ];
X(:,3) = 0;

% Graphical representation
plot(mesh,1)
hold on
plotNrm(mesh,'b')
plot3(X(:,1),X(:,2),X(:,3),'or')
axis equal
view(0,90)
grid on

% Analytical integration (exact)
tic
S     = mesh.vtx;
n     = mesh.nrm;
tau   = mesh.tgt;
[logRa,rlogRa,gradlogRa] = domSemiAnalyticInt2D(X,S,n,tau);
toc

% 1D simpson integration ("exact")
tic
logRs     = zeros(size(logRa));
rlogRs    = zeros(size(rlogRa));
gradlogRs = zeros(size(gradlogRa));
for i = 1:size(X,1)
    % Function |r| in the 2D plane
    fun = @(x1) sqrt( (x1-X(i,1)).^2 + X(i,2).^2 );
    
    % Scalar integration log(|r|)
    logRs(i) = integral(@(x1) log(fun(x1)),a,b,'Reltol',1e-12);
    
    % Vectorial integration rlog(|r|)
    rlogRs(i,1) = integral(@(x1) (x1-X(i,1)).*log(fun(x1)),a,b,'Reltol',1e-12);  
    rlogRs(i,2) = -X(i,2) .*  logRs(i);
    rlogRs(i,3) = 0;
    
    % Vectoral integration gradlog(|r|)
    if (X(i,2) == 0) && (X(i,1)>=a && X(i,1)<=b)
        gradlogRs(i,:) = 0;
    else    
        gradlogRs(i,1) = integral(@(x1) (x1-X(i,1))./fun(x1).^2,a,b);
        gradlogRs(i,2) = integral(@(x1) -X(i,2)./fun(x1).^2,a,b);
        gradlogRs(i,3) = 0;
    end
end
toc

% Relative errors log(|r|)
figure
semilogy( abs(logRa-logRs)./abs(logRs) )
title('Relative error (log) : \int log(|r|)')
grid on
I = 1:N;
disp('Relative error (Linf) : \int log(|r|)')
norm(logRa(I)-logRs(I),'inf')./norm(logRa(I),'inf')
norm(logRa-logRs,'inf')./norm(logRa,'inf')

% Relative errors rlog(|r|)
figure
semilogy( abs(rlogRa-rlogRs)./abs(rlogRs) )
title('Relative error (log) : \int rlog(|r|)')
grid on
I = 1:N;
disp('Relative error (Linf) : \int rlog(|r|)')
norm(rlogRa(I)-rlogRs(I),'inf')./norm(rlogRa(I),'inf')
norm(rlogRa-rlogRs,'inf')./norm(rlogRa,'inf')

% Relative errors gradlog(|r|)
figure
semilogy( abs(gradlogRa-gradlogRs)./abs(gradlogRs) )
title('Relative error (log) : \int gradlog(|r|)')
grid on
I = 1:N;
disp('Relative error (Linf) : \int gradlog(|r|)')
norm(gradlogRa(I)-gradlogRs(I),'inf')./norm(gradlogRa(I),'inf')
norm(gradlogRa-gradlogRs,'inf')./norm(gradlogRa,'inf')




disp('~~> Michto gypsilab !')




