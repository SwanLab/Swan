function M = femUnk2Qud(fe,domain)
%+========================================================================+
%|                                                                        |
%|              OPENFEM - LIBRARY FOR FINITE ELEMENT METHOD               |
%|           openFem is part of the GYPSILAB toolbox for Matlab           |
%|                                                                        |
%| COPYRIGHT : Matthieu Aussal & Francois Alouges (c) 2017-2018.          |
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
%|    #    |   FILE       : femUnknown.m                                  |
%|    #    |   VERSION    : 0.40                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal & François Alouges            |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 14.03.2018                                    |
%| ( === ) |   SYNOPSIS   : Unknowns and reduction matrix for constrained |
%|  `---'  |                finite elements                               |
%+========================================================================+

% Unknown to degrees of freedom matrix 
[~,P] = fe.unk;

% Surfacic restriction (trace)
if (size(fe.msh.elt,2) == size(domain.msh.elt,2) + 1)
    bound  = fe.msh.bnd;
    P      = restriction(fe,bound) * P;
    fe.msh = bound;
end

% Degrees of freedom to quadrature matrix
if strcmp(fe.typ(1),'P')
    M = femLagrangePn(fe,domain);
elseif strcmp(fe.typ,'RWG')
    M = femRaoWiltonGlisson(fe,domain);
elseif strcmp(fe.typ,'NED')
    M = femNedelec(fe,domain);
else
    error('fem.m : unavailable case')
end

% Unknown to quadrature 
if iscell(M)
    M{1} = M{1} * P;
    M{2} = M{2} * P;
    M{3} = M{3} * P;
else
    M = M * P;
end
end
