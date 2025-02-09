function phi = vibsCutoff(n,r,dr)
%+========================================================================+
%|                                                                        |
%|                 OPENVIBS - LIBRARY FOR VIBRO-ACOUSTIC                  |
%|           openVibs is part of the GYPSILAB toolbox for Matlab          |
%|                                                                        |
%| COPYRIGHT : Matthieu Aussal, Marc Bakry (c) 2017-2019.                 |
%| PROPERTY  : Centre de Mathematiques Appliquees, Ecole polytechnique,   |
%| route de Saclay, 91128 Palaiseau, France. All rights reserved.         |
%| LICENCE   : This program is free software, distributed in the hope that|
%| it will be useful, but WITHOUT ANY WARRANTY. Natively, you can use,    |
%| redistribute and/or modify it under the terms of the GNU General Public|
%| License, as published by the Free Software Foundation (version 3 or    |
%| later,  http://www.gnu.org/licenses). For private use, dual licencing  |
%| is available, please contact us to activate a "pay for remove" option. |
%| CONTACT   : matthieu.aussal@polytechnique.edu                          |
%|             marc.bakry@polytechnique.edu                               |
%| WEBSITE   : www.cmap.polytechnique.fr/~aussal/gypsilab                 |
%|                                                                        |
%| Please acknowledge the gypsilab toolbox in programs or publications in |
%| which you use it.                                                      |
%|________________________________________________________________________|
%|   '&`   |                                                              |
%|    #    |   FILE       : vibsCutoff.m                                  |
%|    #    |   VERSION    : 0.55                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal & Marc Bakry                  |
%|  ( # )  |   CREATION   : 14.03.2019                                    |
%|  / 0 \  |   LAST MODIF :                                               |
%| ( === ) |   SYNOPSIS   :                                               |
%|  `---'  |                                                              |
%+========================================================================+
xmin = -r;
xmax = r;
phi  = @(X) (X(:,n) >= xmin) .* (X(:,n) <= xmax) + ...
    (1 + (X(:,n)-xmin)/dr) .* (X(:,n) < xmin) .* (X(:,n) >= xmin-dr) +...
    (1 - (X(:,n)-xmax)/dr) .* (X(:,n) > xmax) .* (X(:,n) <= xmax+dr);
end