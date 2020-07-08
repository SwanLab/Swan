function U = raySource(N)
%+========================================================================+
%|                                                                        |
%|           OPENRAY - LIBRARY FOR TRI-DIMENSIONAL RAY TRACING            |
%|           openRay is part of the GYPSILAB toolbox for Matlab           |
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
%| WEBSITE   : www.cmap.polytechnique.fr/~aussal/gypsilab                 |
%|                                                                        |
%| Please acknowledge the gypsilab toolbox in programs or publications in |
%| which you use it.                                                      |
%|________________________________________________________________________|
%|   '&`   |                                                              |
%|    #    |   FILE       : raySource.m                                   |
%|    #    |   VERSION    : 0.41                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 01.04.2018                                    |
%| ( === ) |   SYNOPSIS   : Spherical quasi-uniform sampling for source   |
%|  `---'  |                                                              |
%+========================================================================+

% Fibonnacci rules for uniform sampling
if N > 1
    or      = (1+sqrt(5))/2;
    theta   = (mod((2*pi/or) .* (0:N-1),2*pi))';
    phi     = asin( -1 + 2/(N-1) * (0:N-1))';
    [x,y,z] = sph2cart(theta,phi,1);
    U       = [x y z];
else
    U = [1 0 0];
end

% Security
if (norm(sqrt(sum(U.^2,2)) - 1,'inf') > 1e-12)
    error('raySource.m : non unitary direction')
end
end
