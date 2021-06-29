function [bool,Xdim,Ydim] = hmxFar(Mh)
%+========================================================================+
%|                                                                        |
%|         OPENHMX - LIBRARY FOR H-MATRIX COMPRESSION AND ALGEBRA         |
%|           openHmx is part of the GYPSILAB toolbox for Matlab           |
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
%|    #    |   FILE       : hmxFar.m                                      |
%|    #    |   VERSION    : 0.54                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 01.05.2019                                    |
%| ( === ) |   SYNOPSIS   : H-Matrix far boolean                          |
%|  `---'  |                                                              |
%+========================================================================+

% Particles box X
X        = Mh.pos{1};
Xmin     = min(X,[],1);
Xmax     = max(X,[],1);
Xctr     = 0.5*(Xmin+Xmax);
Xdgl     = Xmax-Xmin;
[~,Xdim] = max(Xdgl);

% Particles box Y
Y        = Mh.pos{2};
Ymin     = min(Y,[],1);
Ymax     = max(Y,[],1);
Yctr     = 0.5*(Ymin+Ymax);
Ydgl     = Ymax-Ymin;
[~,Ydim] = max(Ydgl);

% Distance admissibility following principal axis (1/2 box)
XYctr = Yctr-Xctr;
XYdgl = Xdgl+Ydgl;
bool  = sum( (abs(XYctr)>=0.75*XYdgl) & (abs(XYctr)>=0.75*mean(XYdgl)) );

% Angle following distance axe regarding from X
U     = [ Y(:,1)-Xctr(:,1) , Y(:,2)-Xctr(:,2) , Y(:,3)-Xctr(:,3)];
V     = Yctr - Xctr;    
alpha = acosd((U*V')./( sqrt(sum(U.^2,2)) .* sqrt(sum(V.^2,2))) );

% Angle following distance axe regarding from Y
U    = [ X(:,1)-Yctr(:,1) , X(:,2)-Yctr(:,2) , X(:,3)-Yctr(:,3)];
V    = Xctr - Yctr;    
beta = acosd((U*V')./( sqrt(sum(U.^2,2)) .* sqrt(sum(V.^2,2))) );

% Add angular admissibility following translation axis (1/2 box)
theta = 30; % acosd(2/sqrt(5))
bool  = bool & (max(alpha)<=theta) & (max(beta)<=theta);
end
