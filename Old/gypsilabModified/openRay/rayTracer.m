function ray = rayTracer(ray,ord,rMax)
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
%|    #    |   FILE       : rayTracer.m                                   |
%|    #    |   VERSION    : 0.41                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 01.04.2018                                    |
%| ( === ) |   SYNOPSIS   : Ray tracer with octree speedup                |
%|  `---'  |                                                              |
%+========================================================================+

% Infos
disp('====> RAY TRACING <====')
tps = time();

% Mesh dimension
Nelt = length(ray.msh);

% Prepare sphere
tree = rayTreeInit(ray);

% Available ray
ind = find( (ray.dst<rMax) & (sum(ray.dir.^2,2)~=0) );

% Current order
n = length(ray.pos)-1;

% Infos
disp([' ~~> Tree with ',num2str(length(tree)),' stage(s) - Elapsed time is ', ...
                num2str(time()-tps),' seconds.'])
            
% Iterative loop
while ~isempty(ind) && (n < ord)
    % Infos
    tps = time();
    
    % Full repartition for smallest meshes
    if (Nelt < 50)   
        Iray = cell(Nelt,1);
        for el = 1:Nelt
            Iray{el} = ind;
        end
        
    % Hierarchical repartition
    else
        Iray = rayTree(ray,tree,ind);
    end
    
    % Intersection between mesh and ray
    ray = rayCollision(ray,Iray);
    
    % Update available ray
    ind = find( (ray.dst<rMax) & (sum(ray.dir.^2,2)~=0) );
    
    % Order incrementation
    n = n + 1;
    
    % Infos
    disp([' + Step ',num2str(n), ' - Elapsed time is ', ...
                num2str(time()-tps),' seconds.'])            
end
end


function tps = time()
tps = clock;
tps = tps(4)*3600 + tps(5)*60 + tps(6);
end


%     % Representation graphique
%     figure
%     hold on
%     for el = 1:Nelt
%         subElt = ray.msh.sub(el);
%         subRay = ray.sub(Iray{el});
%         plot(subElt,1)
%         plot(subRay)
%         pause
%     end
%     hold off
