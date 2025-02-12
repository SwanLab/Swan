function rayPlot(ray,ord)
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
%|    #    |   FILE       : rayPlot.m                                     |
%|    #    |   VERSION    : 0.41                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 01.04.2018                                    |
%| ( === ) |   SYNOPSIS   : Plot trajectory and last direction            |
%|  `---'  |                                                              |
%+========================================================================+

% Graphical configuration
figure(gcf)
off = ~ishold;
if off
    hold on
end

% Plot trajectory
if (ord > -1)
    for i = 2:ord+1
        Pim1 = ray.pos{i-1};
        Pi   = ray.pos{i};
        plot3([Pim1(:,1) Pi(:,1)]',[Pim1(:,2) Pi(:,2)]',[Pim1(:,3) Pi(:,3)]','o-b')
    end
end

% Plot direction
if (ord == -1) || (length(ray.pos) == 1)
    P = ray.pos{end};
    U = ray.dir;
    quiver3(P(:,1),P(:,2),P(:,3),U(:,1),U(:,2),U(:,3),'r')
end

% Hold 
if off
    hold off
end
end
