function [logR,rlogR,gradlogR] = domSemiAnalyticInt2D(X,S,n,tau) 
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
%|    #    |   FILE       : domSemiAnalyticInt2D.m                        |
%|    #    |   VERSION    : 0.53                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal & Martin Averseng             |
%|  ( # )  |   CREATION   : 25.11.2018                                    |
%|  / 0 \  |   LAST MODIF : 14.03.2019                                    |
%| ( === ) |   SYNOPSIS   : Semi-analytic integration on segment for a set|
%|  `---'  |                of particles                                  |
%+========================================================================+

% Dimensions
Nx  = size(X,1);
unx = ones(Nx,1);
un2 = [1 1];

% Parametric coordinates on (AB)
XA = unx*S(1,:)-X;
XB = unx*S(2,:)-X;
d  = -XA*n';      
a  = XA*tau';   
b  = XB*tau';  

% Check distance to (AB)
bool = (abs(d)>1e-8);

% Initialization
logR = zeros(Nx,1);

% General primitive : \int_a^b ln(x^2+d^2) dx
I       = find(bool);
logR(I) = F1(b,d,I) - F1(a,d,I);

% Assymptotic primitive (d=0) : \int_a^b ln(x^2) dx
I       = find(~bool);
logR(I) = F2(b,I) - F2(a,I);

% \int rlogr
rlogRtau = F3(b,d) - F3(a,d);
rlogRn   = -(d.*logR);
rlogR    = rlogRtau*tau + rlogRn*n;

% \int grad(logr)
gradlogRtau       = 1/2*log((b.^2 + d.^2)./(a.^2 + d.^2));
gradlogRtau(I)    = log(abs(b(I)./a(I)));
gradlogRtau(a<=0 & b>=0 & ~bool) = 0;
gradlogRn         = -atan(b./d) + atan(a./d);
gradlogRn(I)      = 0;
gradlogR          = gradlogRtau*tau + gradlogRn*n;
end


function y = F1(x,d,I)
x = x(I);
d = d(I);
y = 1/2*x.*log(x.^2+d.^2) - x + d.*atan(x./d);
end


function y = F2(x,I)
x = x(I);
y = xlog(x) - x;
end


function y = F3(x,d)
y = 1/4 * (xlog(x.^2+d.^2) - x.^2);
end


function y = xlog(x)
y    = x.*log(abs(x));
I    = (abs(x)<1e-15); 
y(I) = 0;
end




