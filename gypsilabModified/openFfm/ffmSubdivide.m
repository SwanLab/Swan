function [chd,Ichd] = ffmSubdivide(X,prt,edg,fig)
%+========================================================================+
%|                                                                        |
%|         OPENFFM - LIBRARY FOR FAST AND FREE MEMORY CONVOLUTION         |
%|           openFfm is part of the GYPSILAB toolbox for Matlab           |
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
%| WEBSITE   : www.cmap.polytechnique.fr/~aussal/gypsilab                 |
%|                                                                        |
%| Please acknowledge the gypsilab toolbox in programs or publications in |
%| which you use it.                                                      |
%|________________________________________________________________________|
%|   '&`   |                                                              |
%|    #    |   FILE       : ffmSubdivide.m                                |
%|    #    |   VERSION    : 0.6                                           |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 05.09.2019                                    |
%| ( === ) |   SYNOPSIS   : Recursive subdivision using an octree         |
%|  `---'  |                                                              |
%+========================================================================+

% Initialisation des enfants par subdivision type Octree
Nprt    = length(prt.ind);
chd.ctr = zeros(8*Nprt,3);
chd.ind = cell(8*Nprt,1);
chd.nbr = zeros(8*Nprt,1);
Ichd    = zeros(Nprt,8);

% Boucle sur les parents
l = 0;
for i = 1:Nprt
    % Subdivision centrale
    ind  = prt.ind{i};
    cutx = X(ind,1)<=prt.ctr(i,1);
    cuty = X(ind,2)<=prt.ctr(i,2);
    cutz = X(ind,3)<=prt.ctr(i,3);
    
    % Repartition des indices
    cut = logical([ ...
        cutx  .*  cuty  .*  cutz , ...
        ~cutx  .*  cuty  .*  cutz , ...
        cutx  .* ~cuty  .*  cutz , ...
        ~cutx  .* ~cuty  .*  cutz , ...
        cutx  .*  cuty  .* ~cutz , ...
        ~cutx  .*  cuty  .* ~cutz , ...
        cutx  .* ~cuty  .* ~cutz , ...
        ~cutx  .* ~cuty  .* ~cutz ] );
    Ncut = sum(cut,1);
    
    % Centres des enfants
    ctr = ones(8,1)*prt.ctr(i,1:3) + 0.25*edg .* [...
        -1 -1 -1 ; 1 -1 -1; -1 1 -1 ; 1 1 -1 ;...
        -1 -1  1 ; 1 -1  1; -1 1  1 ; 1 1  1 ];
    
    % Ajout d'une boite pour chaque enfant non vide
    for j = find(Ncut)
        chd.ctr(l+1,:) = ctr(j,:);
        chd.ind{l+1}   = ind(cut(:,j));
        chd.nbr(l+1)   = Ncut(j);
        Ichd(i,j)      = l+1;
        l = l + 1;
    end
end

% Extraction des boites vides
chd.ctr = chd.ctr(1:l,:);
chd.ind = chd.ind(1:l);
chd.nbr = chd.nbr(1:l);

% Verification que chaque pot a son couvercle
if abs( sum(chd.nbr) - sum(prt.nbr) ) > 1e-12
    error('ffmSubDivide.m - error 1')
end
if abs( sum(cell2mat(chd.ind)) - sum(cell2mat(prt.ind)) ) > 1e-12
    error('ffmSubDivide.m - error 2')
end

% Representation graphique
if fig
    % Configuration
    figure(fig); clf
    hold on
    title('Box plot')
    grid on ; axis equal
    
    % Representation des centres
    plot3(chd.ctr(:,1),chd.ctr(:,2),chd.ctr(:,3),'ok','MarkerSize',8)
    
    % Representation des points par boites
    for i = 1:length(chd.ind)
        plot3(X(chd.ind{i},1),X(chd.ind{i},2),X(chd.ind{i},3),'.','Color',rand(1,3))
    end
end
end
