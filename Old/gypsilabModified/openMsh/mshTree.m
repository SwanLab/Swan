function tree = mshTree(mesh,typ,Nlf,fig)
%+========================================================================+
%|                                                                        |
%|                 OPENMSH - LIBRARY FOR MESH MANAGEMENT                  |
%|           openMsh is part of the GYPSILAB toolbox for Matlab           |
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
%|    #    |   FILE       : mshtree.m                                     |
%|    #    |   VERSION    : 0.41                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 01.04.2018                                    |
%| ( === ) |   SYNOPSIS   : Hierarchical tree subdivision                 |
%|  `---'  |                                                              |
%+========================================================================+

% Security
if isempty(Nlf)
    Nlf = 1;
end

% Particles are element centers
X  = mesh.ctr;
Nx = length(mesh);

% Box initialization
Xmin = min(X,[],1);
Xmax = max(X,[],1);
edg  = max(Xmax-Xmin);

% Initialize octree
tree           = cell(1,1);
tree{1}.chd    = [];
tree{1}.ctr    = Xmin + 0.5*edg;
tree{1}.edg    = edg;
tree{1}.ind{1} = (1:Nx)';
tree{1}.nbr    = Nx;
tree{1}.prt    = 0;

% Recursive hierarchical clustering
n = 1;
while (max(tree{n}.nbr) > Nlf)
    % Tree subdivision
    if strcmp(typ,'octree')
        [tree{n+1,1},Ichd] = mshSubDivide8(X,tree{n});
    elseif strcmp(typ,'binary')
        [tree{n+1,1},Ichd] = mshSubDivide2(X,tree{n});
    else
        error('mshTree.m : unavailable case')
    end
        
    % Add children indices to parents
    tree{n}.chd = Ichd; 

    % Verify children
    if (norm( cell2mat(tree{n}.chd) - (1:length(tree{n+1}.nbr))' , 'inf') > 1e-12)
        error('mshTree.m : unavailable case')
    end
    
    % Verify indices
    if (norm( sort(cell2mat(tree{n+1}.ind)) - sort(cell2mat(tree{n+1}.ind)) , 'inf') > 1e-12)
        error('mshTree.m : unavailable case')
    end
    
    % Verify number
    if (abs( sum(tree{n+1}.nbr) - sum(tree{n}.nbr) ) > 1e-12)
        error('mshTree.m : unavailable case')
    end
    
    % Verify parents
    if (norm( unique(tree{n+1}.prt) - (1:length(tree{n}.nbr))' , 'inf') > 1e-12)
        error('mshTree.m - unavailable case')
    end
    
    % Incrementation
    n = n+1;
    
    % Graphical representation
    if fig
        % Numbering boxes
        V = zeros(Nx,1);
        for i = 1:length(tree{n}.ind)
            V(tree{n}.ind{i}) = i;
        end
        
        % Graphical representation
        figure(fig); clf
        plot3(tree{n}.ctr(:,1),tree{n}.ctr(:,2),tree{n}.ctr(:,3),'ok','MarkerSize',8)
        hold on
        plot(mesh,V) 
        title(['Tree at step ',num2str(n-1)])
        grid on ; axis equal
        xlabel('X'); zlabel('Y'); zlabel('Z');
        alpha(0.9)
        colorbar
        drawnow
    end
end
end


function [chd,Ichd] = mshSubDivide8(X,prt)
% Initialize children tree 
Nprt    = length(prt.ind);
chd.chd = [];
chd.ctr = zeros(8*Nprt,3);
chd.edg = 0.5 * prt.edg;
chd.ind = cell(8*Nprt,1);
chd.nbr = zeros(8*Nprt,1);
chd.prt = zeros(8*Nprt,1);

% Initialize children indices
Ichd = cell(Nprt,1);

% Loop for parents
l = 0;
for i = 1:Nprt
    % Central subdivision
    ind  = prt.ind{i};
    cutx = X(ind,1)<=prt.ctr(i,1);
    cuty = X(ind,2)<=prt.ctr(i,2);
    cutz = X(ind,3)<=prt.ctr(i,3);
    
    % Indices repartition for children
    cut = [ ...
        cutx   &  cuty  &  cutz , ...
        ~cutx  &  cuty  &  cutz , ...
        cutx   & ~cuty  &  cutz , ...
        ~cutx  & ~cuty  &  cutz , ...
        cutx   &  cuty  & ~cutz , ...
        ~cutx  &  cuty  & ~cutz , ...
        cutx   & ~cuty  & ~cutz , ...
        ~cutx  & ~cuty  & ~cutz ];
    Ncut = sum(cut,1);
    
    % Indices for parents
    Ichd{i} = l + (1:sum(Ncut>0))';
    
    % Children centers
    ctr = ones(8,1)*prt.ctr(i,:) + 0.25*prt.edg .* [...
        -1 -1 -1 ; 1 -1 -1; -1 1 -1 ; 1 1 -1 ;...
        -1 -1  1 ; 1 -1  1; -1 1  1 ; 1 1  1 ];
    
    % Add box only for non empty children
    for j = find(Ncut) 
        % Children data
        chd.ctr(l+1,:) = ctr(j,:);
        chd.ind{l+1}   = ind(cut(:,j));
        chd.nbr(l+1)   = Ncut(j);
        chd.prt(l+1)   = i;
        
        % Incrementation
        l = l + 1;
    end
end

% Extract non empty box
chd.ctr = chd.ctr(1:l,:);
chd.ind = chd.ind(1:l);
chd.nbr = chd.nbr(1:l);
chd.prt = chd.prt(1:l);
end


function [chd,Ichd] = mshSubDivide2(X,prt)
% Initialize children tree 
Nprt    = length(prt.ind);
chd.chd = [];
chd.ctr = zeros(2*Nprt,3);
chd.edg = zeros(2*Nprt,1);
chd.ind = cell(2*Nprt,1);
chd.nbr = zeros(2*Nprt,1);
chd.prt = zeros(2*Nprt,1);

% Initialize children indices
Ichd = cell(Nprt,1);

% Loop for parents
l = 0;
for i = 1:Nprt
    % Local points
    ind = prt.ind{i};
    
    % Local box
    xmin = min(X(ind,:),[],1);
    xmax = max(X(ind,:),[],1);
    xdgl = xmax - xmin;
    
    % Median subdivision over largest dimension
    [~,d] = max(xdgl);
    m     = median(X(ind,d));
    cutd  = (X(ind,d) <= m);
    cut   = [ cutd , ~cutd ];
    Ncut  = sum(cut,1);
    
    % Security for planar repartition
    if (abs(Ncut(1)-Ncut(2)) > 1)
        [~,I]   = sort(X(ind,d));
        cutd(:) = 0;
        cutd(I(1:floor(end/2))) = 1;
        cut  = [ cutd , ~cutd ];
        Ncut = sum(cut,1);
    end
    
    % Indices for parents
    Ichd{i} = l + (1:sum(Ncut>0))';
    
    % Add box only for non empty children
    for j = find(Ncut) 
        % Local box
        Ij   = ind(cut(:,j));
        xmin = min(X(Ij,:),[],1);
        xmax = max(X(Ij,:),[],1);

        % Children data
        chd.ctr(l+1,:) = 0.5*(xmin+xmax);
        chd.edg(l+1)   = max(xmax-xmin);
        chd.ind{l+1}   = Ij;
        chd.nbr(l+1)   = Ncut(j);
        chd.prt(l+1)   = i;
        
        % Incrementation
        l = l + 1;
    end
end

% Extract non empty box
chd.ctr = chd.ctr(1:l,:);
chd.edg = chd.edg(1:l);
chd.ind = chd.ind(1:l);
chd.nbr = chd.nbr(1:l);
chd.prt = chd.prt(1:l);
end
