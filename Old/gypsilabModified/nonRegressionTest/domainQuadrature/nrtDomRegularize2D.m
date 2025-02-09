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
%|    #    |   FILE       : nrtDomRegularize2D.m                          |
%|    #    |   VERSION    : 0.50                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 25.11.2018                                    |
%|  / 0 \  |   LAST MODIF : 25.11.2018                                    |
%| ( === ) |   SYNOPSIS   : Singular kernel regularization (in debug mode)|
%|  `---'  |                                                              |
%+========================================================================+

% Cleaning
clear all
close all
clc

% Gypsilab path
run('../../addpathGypsilab.m')

% Parameters
N   = 1e2;
gss = 3;

% Spherical mesh
circle = mshCircle(N,1);

% Square mesh
square = mshSquare(2*N,[3 3]);

% Graphical representation
figure
plot(circle)
hold on
plot(square,'w')
axis equal
alpha(0.5)

% Domain
sigma = dom(circle,gss);    

% Finite elements
u = fem(circle,'P0');
v = fem(circle,'P1');

% Gren kernel
Gxy      = @(X,Y) femGreenKernel(X,Y,'[log(r)]',[]);
dyGxy{1} = @(X,Y) femGreenKernel(X,Y,'grady[log(r)]1',[]);
dyGxy{2} = @(X,Y) femGreenKernel(X,Y,'grady[log(r)]2',[]);
dyGxy{3} = @(X,Y) femGreenKernel(X,Y,'grady[log(r)]3',[]);

% Single radiation P0
M  = integral(square.vtx,sigma,Gxy,u);
Mr = regularize(square.vtx,sigma,'[log(r)]',u);
norm(M+Mr,'inf')

% Double radiation P0
M  = integral(square.vtx,sigma,dyGxy,ntimes(u));
Mr = regularize(square.vtx,sigma,'grady[log(r)]',ntimes(u));
norm(M+Mr,'inf')

% Single radiation P1
M  = integral(square.vtx,sigma,Gxy,v);
Mr = regularize(square.vtx,sigma,'[log(r)]',v);
norm(M+Mr,'inf')

% Single layer P0
M  = integral(sigma,sigma,u,Gxy,u);
Mr = regularize(sigma,sigma,u,'[log(r)]',u);
norm(M+Mr,'inf')

% Double layer P0
M  = integral(sigma,sigma,u,dyGxy,ntimes(u));
Mr = regularize(sigma,sigma,u,'grady[log(r)]',ntimes(u));
norm(M+Mr,'inf')

% Single layer P1
M  = integral(sigma,sigma,u,Gxy,v);
Mr = regularize(sigma,sigma,u,'[log(r)]',v);
norm(M+Mr,'inf')

% Hypersingular P1
M  = integral(sigma,sigma,ntimes(u),Gxy,ntimes(v));
Mr = regularize(sigma,sigma,ntimes(u),'[log(r)]',ntimes(v));
norm(M+Mr,'inf')

M  = integral(sigma,sigma,nxgrad(u),Gxy,nxgrad(v));
Mr = regularize(sigma,sigma,nxgrad(u),'[log(r)]',nxgrad(v));
norm(M+Mr,'inf')



disp('~~> Michto gypsilab !')




