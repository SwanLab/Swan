function err = domDifference(domain,fe,Uh,Uex,type )
%+========================================================================+
%|                                                                        |
%|              OPENDOM - LIBRARY FOR NUMERICAL INTEGRATION               |
%|           openDom is part of the GYPSILAB toolbox for Matlab           |
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
%|    #    |   FILE       : domDifference.m                               |
%|    #    |   VERSION    : 0.40                                          |
%|   _#_   |   AUTHOR(S)  : François Alouges                              |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 14.03.2018                                    |
%| ( === ) |   SYNOPSIS   : L2 and H1 ERRORS                              |
%|  `---'  |                                                              |
%+========================================================================+

% Mesh integration data
[X,Wx] = domain.qud;

% Finite element matrix
Mu = fe.uqm(domain);

% Function to be applied
Uexact = Uex(X);
if (size(Uexact,2) > 1)
    error('domDifference.m : unavailable case')
end
    
Uapp = Mu * Uh;
switch type
    case 'L2'
        err = sqrt(sum(Wx.*(Uexact-Uapp).^2));
    case 'H1'
        eps = 1e-6;
        gradef = grad(fe);
        Gu = gradef.uqm(domain);
        DxUapp = Gu{1} * Uh;
        DyUapp = Gu{2} * Uh;
        DzUapp = Gu{3} * Uh;
        n1 = size(X,1);
        eps1 = eps*ones(n1,1)*[1,0,0];
        eps2 = eps*ones(n1,1)*[0,1,0];
        eps3 = eps*ones(n1,1)*[0,0,1];
        DxUexact = (Uex(X+eps1)-Uex(X-eps1))/(2*eps);
        DyUexact = (Uex(X+eps2)-Uex(X-eps2))/(2*eps);
        DzUexact = (Uex(X+eps3)-Uex(X-eps3))/(2*eps);
        err = sqrt(sum(Wx.*((Uexact-Uapp).^2+(DxUexact-DxUapp).^2+(DyUexact-DyUapp).^2+(DzUexact-DzUapp).^2)));
    otherwise
        error('Unknown error type. Known types are ''L2'' or ''H1''.');
end
end
