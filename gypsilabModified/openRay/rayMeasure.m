function [I,src] = rayMeasure(ray,Xmes,rad,rMax)
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
%|    #    |   FILE       : rayMeasure.m                                  |
%|    #    |   VERSION    : 0.41                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 01.04.2018                                    |
%| ( === ) |   SYNOPSIS   : Spherical measurement                         |
%|  `---'  |                                                              |
%+========================================================================+

% Initialization
I   = cell(1,length(ray.pos)-1);
src = cell(length(ray.pos)-1,1);
dst = zeros(length(ray),1);

% Loop 
for i = 1:size(I,2)
    % Ray positions
    Pi   = ray.pos{i};
    Pip1 = ray.pos{i+1};
    
    % Ray direction and length
    U   = (Pip1 - Pi);
    lgt = sqrt(sum(U.^2,2));
    U   = U ./ (lgt * [1 1 1]);

    % Initial position -> micro
    PM = ones(size(Pi,1),1)*Xmes - Pi;
    
    % Measure triangle (pythagore)
    hyp2 = sum(PM.^2,2);
    adj  = sum(PM.*U,2);
    opp2 = hyp2 - adj.^2;
    
    % Measured ray
    I{i} = find( (hyp2 < lgt.^2) & (adj >= 0) & (opp2 <= rad^2) & (dst + sqrt(hyp2) < rMax) );

    % Sources (focusing)
    if ~isempty(I{i})
        src{i} = Pi(I{i},:) - (dst(I{i})*[1 1 1]) .* U(I{i},:);
    else
        src{i} = zeros(0,3);
    end
    
    % Total distance from initial position
    dst = dst + lgt;    
end
end
