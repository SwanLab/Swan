function [M,J] = oprSubdivideCell(M,I,type)
%+========================================================================+
%|                                                                        |
%|            OPENOPR - LIBRARY FOR SPECIFIC OPERATORS IN BEM             |
%|           openOpr is part of the GYPSILAB toolbox for Matlab           |
%|                                                                        |
%| COPYRIGHT : Matthieu Aussal (c) 2017-2019.                             |
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
%|    #    |   FILE       : oprSubdivideCell.m                            |
%|    #    |   VERSION    : 0.61                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 05.09.2019                                    |
%| ( === ) |   SYNOPSIS   : Subdivide finite elements cells and matrix    |
%|  `---'  |                                                              |
%+========================================================================+

% Left integration
if strcmp(type,'left')
    % Left integration indices
    if iscell(M)
        Vx = 0;
        for n = 1:length(M)
           Vx = Vx + (1+rand(1,length(I))) * M{n}(I,:); 
        end
    else
        Vx = (1+rand(1,length(I))) * M(I,:);
    end
    J = find(Vx);
    if isempty(J)
        J = [1;2]';
    end
    
    % Left matrix subdivision
    if iscell(M)
        for n = 1:length(M)
            M{n} = M{n}(I,J);
        end
    else
        M = M(I,J);
    end

% Right integration
elseif strcmp(type,'right')
    % Right integration indices
    if iscell(M)
        Vy = 0;
        for n = 1:length(M)
            Vy = Vy + M{n}(:,I) * (1+rand(length(I),1));
        end
    else
        Vy = M(:,I) * (1+rand(length(I),1));
    end
    J = find(Vy);
    if isempty(J)
        J = [1;2]';
    end
    
    % Right matrix subdivision
    if iscell(M)
        for n = 1:length(M)
            M{n} = M{n}(J,I);
        end
    else
        M = M(J,I);
    end

% Unknown case     
else
    error('femSubdivideCell.m : unavailable case');
end
end
