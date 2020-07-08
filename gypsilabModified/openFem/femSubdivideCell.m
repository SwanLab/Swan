function [M,J] = femSubdivideCell(M,I,type)
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
%|    #    |   FILE       : femSubdivideCell.m                            |
%|    #    |   VERSION    : 0.40                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal & François Alouges            |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 14.03.2018                                    |
%| ( === ) |   SYNOPSIS   : Subdivide finite elements cells and matrix    |
%|  `---'  |                                                              |
%+========================================================================+

% Left integration
if strcmp(type,'left')
    % Left integration indices
    if iscell(M)
        Vx = (1+rand(1,length(I))) * M{1}(I,:) + ...
            (1+rand(1,length(I))) * M{2}(I,:) + ...
            (1+rand(1,length(I))) * M{3}(I,:);
    else
        Vx = (1+rand(1,length(I))) * M(I,:);
    end
    J = find(Vx);
    if isempty(J)
        J = [1;2]';
    end
    
    % Left matrix subdivision
    if iscell(M)
        M{1} = M{1}(I,J);
        M{2} = M{2}(I,J);
        M{3} = M{3}(I,J);
    else
        M = M(I,J);
    end

% Right integration
elseif strcmp(type,'right')
    % Right integration indices
    if iscell(M)
        Vy = M{1}(:,I) * (1+rand(length(I),1)) + ...
            M{2}(:,I) * (1+rand(length(I),1)) + ...
            M{3}(:,I) * (1+rand(length(I),1)) ;
    else
        Vy = M(:,I) * (1+rand(length(I),1));
    end
    J = find(Vy);
    if isempty(J)
        J = [1;2]';
    end
    
    % Right matrix subdivision
    if iscell(M)
        M{1} = M{1}(J,I);
        M{2} = M{2}(J,I);
        M{3} = M{3}(J,I);
    else
        M = M(J,I);
    end

% Unknown case     
else
    error('femSubdivideCell.m : unavailable case');
end
end
