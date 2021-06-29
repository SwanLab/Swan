function Iray = rayTree(ray,tree,Iray)
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
%|    #    |   FILE       : rayTree.m                                     |
%|    #    |   VERSION    : 0.41                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 01.04.2018                                    |
%| ( === ) |   SYNOPSIS   : Ray propagation along mesh tree               |
%|  `---'  |                                                              |
%+========================================================================+

% Intialize parent ray indices
Iprt{1} = Iray;

% Loop at each step of the tree
for n = 2:length(tree)
    % Initialize children ray indices
    Ichd = cell(size(tree{n}.ind));
    
    % Children intersection
    for i = 1:length(Ichd)
        % Parent ray indices
        ir = Iprt{tree{n}.prt(i)}; 
        
        % For non-empty indices
        if ~isempty(ir)
            % Measurement sphere
            M    = tree{n}.ctr(i,:);
            rad2 = tree{n}.rad(i)^2;
            
            % Positions initiales et finales des rayons
            PM = ones(length(ir),1)*M - ray.pos{end}(ir,:);
            U  = ray.dir(ir,:);
            
            % Pythagore measurement
            hyp2 = sum(PM.^2,2);
            adj  = sum(PM.*U,2);
            opp2 = hyp2 - adj.^2;
            
            % Measured ray
            Ichd{i} = ir( (hyp2 < rad2) | ((adj >= 0) & (opp2 <= rad2)) );
        end
    end
    
    % Incrementation
    Iprt = Ichd;
end

% Sort ray indices as elements
Iray       = Iprt;
Ielt       = cell2mat(tree{end}.ind);
Iray(Ielt) = Iray;
end
