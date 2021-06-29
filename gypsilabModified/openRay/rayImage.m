function [img,nrg] = rayImage(ray,mic,rad,rMax)
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
%|    #    |   FILE       : rayImage.m                                    |
%|    #    |   VERSION    : 0.41                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 01.04.2018                                    |
%| ( === ) |   SYNOPSIS   : Image-sources from ray-tracing                |
%|  `---'  |                                                              |
%+========================================================================+

% Material properties (http://www.odeon.dk/material-manufactures)
load('odeon.mat')
mat  = odeon.mat(ray.msh.col,:);
Nfrq = size(mat,2);

% Air absorbing formulae for T = 20°C (bouquin Jouhaneau p 68-69)
rhm = 30;
air = 5.5 * (50/rhm) .* (odeon.frq/1000).^1.7 * 1e-4;

% Measures
[Isrc,src] = ray.measure(mic,rad,rMax);

% Initialization
img = cell(size(src));
nrg = cell(size(src));
rfl = ones(length(ray),Nfrq);

% Loop on images sources
for i = 1:length(src)
    % Unicity for images
    [~,Ia,Ic] = unique(round(src{i}*1e12),'rows','stable');
    
    % Image coordinates
    img{i} = src{i}(Ia,:) - ones(length(Ia),1)*mic;
    
    % Image Energy
    nrg{i} = zeros(length(Ia),Nfrq);
    for j = 1:Nfrq
        nrg{i}(:,j) = accumarray(Ic,rfl(Isrc{i},j),[length(Ia),1]);       
    end
    
    % Update reflecting coeff
    ind         = (ray.iel(:,i+1) > 0);
    rfl(ind,:)  = rfl(ind,:) .* (1 - mat(ray.iel(ind,i+1),:));    
    rfl(~ind,:) = 0;
end

% Vectoriel format
img = cell2mat(img);
nrg = cell2mat(nrg);

% Atmosphere dissipation
dst = sqrt(sum(img.^2,2));
nrg = nrg .* exp(-dst*air);

% Sort in phase
[~,ind] = sort(dst);
img     = img(ind,:);
nrg     = nrg(ind,:);

% Normalisation
nrg = nrg ./ max(max(nrg));      
end
