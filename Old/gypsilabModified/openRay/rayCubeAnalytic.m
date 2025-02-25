function [img,nrg] = rayCubeAnalytic(L,mat,air,Xsrc,Xmic,Rmax)
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
%|    #    |   FILE       : rayCubeAnalytic.m                             |
%|    #    |   VERSION    : 0.41                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 01.04.2018                                    |
%| ( === ) |   SYNOPSIS   : Analytical source image method for a cube     |
%|  `---'  |                                                              |
%+========================================================================+

% Indices
i = -30:30;            
a = i + 0.5 - 0.5*(-1).^i;
b = (-1).^i;

% Position relatives des sources images
x       = b*Xsrc(1) + a*L(1) - Xmic(1);
y       = b*Xsrc(2) + a*L(2) - Xmic(2);
z       = b*Xsrc(3) + a*L(3) - Xmic(3);
[x,y,z] = meshgrid(x,y,z); 
img     = [x(:),y(:),z(:)];

% Dissipation de l'energie par propagation spherique
dst = sqrt(sum(img.^2,2));
nrg = 1./(dst.^2);

% Dissipation de l'energie par les parois
nrg   = nrg * ones(1,size(mat,2));
rfl   = (1-mat);
i0    = abs(0.5*i - 0.25 + 0.25*(-1).^i);
i1    = abs(0.5*i + 0.25 - 0.25*(-1).^i);
for j = 1:size(nrg,2)
    [rx,ry,rz] = meshgrid(...
        (rfl(1,j).^i0).*(rfl(2,j).^i1),...
        (rfl(3,j).^i0).*(rfl(4,j).^i1),...
        (rfl(5,j).^i0).*(rfl(6,j).^i1));
    nrg(:,j) = (rx(:).*ry(:).*rz(:)) .* nrg(:,j);
end

% Energie dissipee par l'air en fonction de la distance
dst = sqrt(sum(img.^2,2));
nrg = nrg .* exp(-dst*air);

% Selections des images selon l'ordre initial
ind = find(dst<Rmax);
dst = dst(ind);
img = img(ind,:);
nrg = nrg(ind,:);

% Sort in phase
[~,ind] = sort(dst);
img       = img(ind,:);
nrg       = nrg(ind,:);

% Normalisation de l'energie
nrg = nrg./max(max(nrg));

end
