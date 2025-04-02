%+========================================================================+
%|                                                                        |
%|            This script uses the GYPSILAB toolbox for Matlab            |
%|                                                                        |
%| COPYRIGHT : Matthieu Aussal (c) 2015-2017.                             |
%| PROPERTY  : Centre de Mathematiques Appliquees, Ecole polytechnique,   |
%| route de Saclay, 91128 Palaiseau, France. All rights reserved.         |
%| LICENCE   : This program is free software, distributed in the hope that|
%| it will be useful, but WITHOUT ANY WARRANTY. Natively, you can use,    |
%| redistribute and/or modify it under the terms of the GNU General Public|
%| License, as published by the Free Software Foundation (version 3 or    |
%| later,  http://www.gnu.org/licenses). For private use, dual licencing  |
%| is available, please contact us to activate a "pay for remove" option. |
%| CONTACT   : matthieu.aussal@polytechnique.edu                          |
%| WEBSITE   : www.cmap.polytechnique.fr/~aussal/gypsilab                 |
%|                                                                        |
%| Please acknowledge the gypsilab toolbox in programs or publications in |
%| which you use it.                                                      |
%|________________________________________________________________________|
%|   '&`   |                                                              |
%|    #    |   FILE       : nrtFemAllWithShuffle.m                        |
%|    #    |   VERSION    : 0.42                                          |
%|   _#_   |   AUTHOR(S)  : Francois Alouges                              |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 15.07.2018                                    |
%| ( === ) |   SYNOPSIS   : Evaluate all finite elements matrices         |
%|  `---'  |                                                              |
%+========================================================================+

% Nettoyage
clear all
close all
clc

% Gypsilab path
run('../../addpathGypsilab.m')

% Create a 2D and a 3D mesh
m1 = mshCube(100,[1 1 1]);
m2 = m1.bnd;
omega1 = dom(m1,4);
omega2 = dom(m2,3);

% Shuffle the meshes
m1b = shuffle(m1);
m1c = shuffle(m1);

%%%%% P0 3D %%%%%
ERREUR = compare(omega1,m1b,m1c,'P0')

%%%%% P1 3D %%%%%
ERREUR = compare(omega1,m1b,m1c,'P1')

%%%%% P2 3D %%%%%
ERREUR = compare(omega1,m1b,m1c,'P2')

m1b = shuffle(m1,'elt');
m2b = shuffle(m2,'elt');
m1c = shuffle(m1,'elt');
m2c = shuffle(m2,'elt');

%%%%% RWG 3D %%%%%
ERREUR = compare(omega1,m1b,m1c,'RWG')

%%%%% RWG 2D %%%%%
ERREUR = compare(omega2,m2b,m2c,'RWG')

%%%%% NED 3D %%%%%
ERREUR = compare(omega1,m1b,m1c,'NED')

%%%%% NED 2D %%%%%
ERREUR = compare(omega2,m2b,m2c,'NED')

%%%%% nx(NED) = RWG 2D %%%%%
Vh = fem(m2b,'NED');
Vhb = fem(m2c,'RWG');
N1  = size(Vh.unk,1);
N2  = size(Vhb.unk,1);

% Compare the matrices
[~,I1,I2] = intersect(Vh.unk,Vhb.unk,'rows');
P        = sparse(I1,I2,1,N1,N2);
f = @(Y) Y(:,2);
M1 = integral(omega2,f,nx(Vh));
M2 = integral(omega2,f,Vhb);
ERREUR = max([norm(M1{1}-M2{1}*P),norm(M1{2}-M2{2}*P),...
    norm(M1{3}-M2{3}*P)])

%%%%% nx(NED(3D)) = RWG 2D %%%%%
Vh = fem(m2c,'RWG');
Vhb = fem(m1b,'NED');
N1  = size(Vh.unk,1);
N2  = size(Vhb.unk,1);
% Compare the matrices
[~,I1,I2] = intersect(Vh.unk,Vhb.unk,'rows');
P        = sparse(I1,I2,1,N1,N2);
f = @(Y) Y(:,2);
M1 = integral(omega2,f,nx(Vhb));
M2 = integral(omega2,f,Vh);
ERREUR = max([norm(M1{1}-M2{1}*P),norm(M1{2}-M2{2}*P),...
    norm(M1{3}-M2{3}*P)])




disp('~~> Michto gypsilab !')



