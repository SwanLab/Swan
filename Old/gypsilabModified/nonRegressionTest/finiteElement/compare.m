function err = compare(omega,mesh1,mesh2,typ)
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
%|    #    |   FILE       : compare.m                                     |
%|    #    |   VERSION    : 0.42                                          |
%|   _#_   |   AUTHOR(S)  : Francois Alouges                              |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 15.07.2018                                    |
%| ( === ) |   SYNOPSIS   : Compare mesh and finite elements matrices     |
%|  `---'  |                                                              |
%+========================================================================+

% Finite element
Vh  = fem(mesh1,typ);
Vhb = fem(mesh2,typ);
N1  = size(Vh.unk,1);

% Compare the mass matrices
[~,I1,I2] = intersect(Vh.unk,Vhb.unk,'rows');
P        = sparse(I1,I2,1,N1,N1);

% Mass matrices
Mass1 = integral(omega,Vh,Vh);
Mass2 = integral(omega,Vhb,Vhb);

% Error Linf
err = max(max(abs(Mass2 - P'*Mass1*P)));
end