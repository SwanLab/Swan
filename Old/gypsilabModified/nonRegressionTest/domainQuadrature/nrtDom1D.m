%+========================================================================+
%|                                                                        |
%|            This script uses the GYPSILAB toolbox for Matlab            |
%|                                                                        |
%| COPYRIGHT : Francois Alouges (c) 2017-2018.                            |
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
%|    #    |   FILE       : nrtDom1D.m                                    |
%|    #    |   VERSION    : 0.40                                          |
%|   _#_   |   AUTHOR(S)  : Francois Alouges                              |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 14.03.2018                                    |
%| ( === ) |   SYNOPSIS   : Polynomial integration                        |
%|  `---'  |                                                              |
%+========================================================================+

% Cleaning
clear all
close all
clc

% Gypsilab path
run('../../addpathGypsilab.m')

%for nbPts=10:3:50
nbPts=10;

vtx = [0 0 0; 1 0 0];
elt = [1 2];

mesh  = msh(vtx, elt);
omega = dom(mesh);

f0=@(x) ones(size(x,1),1);
f1=@(x) x(:,1);
f2=@(x) x(:,1).^2;
f3=@(x) x(:,1).^3;
f4=@(x) x(:,1).^4;
f5=@(x) x(:,1).^5;
f6=@(x) x(:,1).^6;
f7=@(x) x(:,1).^7;
f8=@(x) x(:,1).^8;
f9=@(x) x(:,1).^9;

for i=1:5
    omega.gss = i;
    I0 = integral(omega, f0) - 1 
    I1 = 2*integral(omega, f1) - 1 
    I2 = 3*integral(omega, f2) - 1 
    I3 = 4*integral(omega, f3) - 1 
    I4 = 5*integral(omega, f4) - 1 
    I5 = 6*integral(omega, f5) - 1 
    I6 = 7*integral(omega, f6) - 1 
    I7 = 8*integral(omega, f7) - 1 
    I8 = 9*integral(omega, f8) - 1 
    I9 = 10*integral(omega, f9) - 1 
end

disp('~~> Michto gypsilab !')
