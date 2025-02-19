function [meshr,Ir] = mshMidpoint(mesh,I)
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
%|    #    |   FILE       : mshMidpoint.m                                 |
%|    #    |   VERSION    : 0.50                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 25.11.2018                                    |
%| ( === ) |   SYNOPSIS   : Refine triangular mesh with midpoint algorithm|
%|  `---'  |                                                              |
%+========================================================================+

% Check dimenion
if (size(mesh,2) ~= 3)
    error('mshMidpoint : unavailable case 1')
end

% Save color and replace by hierarchy
col      = mesh.col;
mesh.col = (1:length(mesh))';

% Submeshing (triangle)
meshs = mesh.sub(I);

% Edge meshes
edgs       = meshs.edg;
[edg,Itri] = mesh.edg;

% Interface with edge multiplicity for triangle
[int,Iedg] = intersect(edg,edgs);
ind        = ismember(Itri,Iedg);
mlt        = sum(ind,2);

% Security
tmp = setdiff(int,edgs);
if (size(tmp.elt,1) ~= 0)
    error('mshMidpoint.m : unavailable case 2');
end

% Initialize refined mesh for element without refinement
meshr = mesh.sub(mlt==0);

% Subdivision with 1 common edge
tmp   = mesh.sub(mlt==1);
tmp   = mshMidpoint1(tmp,int);
meshr = union(meshr,tmp);

% Subdivision with 2 common edges
tmp   = mesh.sub(mlt==2);
tmp   = mshMidpoint2(tmp,int);
meshr = union(meshr,tmp);

% Subdivision with 3 common edges
tmp   = mesh.sub(mlt==3);
tmp   = mshMidpoint3(tmp);
meshr = union(meshr,tmp);

% Parent indices and replace colours
Ir        = meshr.col;
meshr.col = col(Ir);

% Security
if (sum(mesh.ndv)-sum(meshr.ndv))/sum(mesh.ndv) > 1e-15*length(meshr)
    error('mshMidpoint.m : unavailable case 3');
end
end


function mesh = mshMidpoint1(mesh,int)
% Mesh nodes and edges center
[nds,ctr] = data(mesh);

% Interface center
Xctr = int.ctr;

% Refined mesh initialization
Nvtx = size(mesh.vtx,1);
Nelt = length(mesh);
col  = mesh.col;
mesh = mesh.sub([]);

% Loop over nodes
for i = 1:3
    % Neighbours
    ip1 = mod(i,3)+1;
    ip2 = mod(ip1,3)+1;   

    % Selected center are inside subdivided mesh 
    I = find(ismember(single(ctr{i}),single(Xctr),'rows'));
    
    % First elements
    vtx  = [nds{i}(I,:) ; nds{ip1}(I,:) ; ctr{i}(I,:)];
    elt  = reshape((1:3*length(I))',length(I),3);
    tmp  = msh(vtx,elt,col(I));
    mesh = union(mesh,tmp);
    
    % Second elements
    vtx  = [nds{i}(I,:) ; ctr{i}(I,:) ; nds{ip2}(I,:)];
    elt  = reshape((1:3*length(I))',length(I),3);
    tmp  = msh(vtx,elt,col(I));
    mesh = union(mesh,tmp);

    % Mesh fusion with previous submeshes
    mesh = union(mesh,tmp);
end

% Security
if size(mesh.elt,1) ~= 2*Nelt
    error('mshMidpoint1.m : unavailable case 1')
end
if size(mesh.vtx,1) ~= Nvtx+Nelt
    error('mshMidpoint1.m : unavailable case 2')
end
end


function mesh = mshMidpoint2(mesh,int)
% Mesh nodes and edges center
[nds,ctr] = data(mesh);

% Interface center
Xctr = int.ctr;

% Refined mesh initialization
Nvtx = size(mesh.vtx,1);
Nelt = length(mesh);
col  = mesh.col;
mesh = mesh.sub([]);

% Loop over nodes
for i = 1:3
    % Neighbours
    ip1 = mod(i,3)+1;
    ip2 = mod(ip1,3)+1;   

    % Selected nodes center not inside subdivided mesh 
    I = find(~ismember(single(ctr{i}),single(Xctr),'rows'));
    
    % First elements
    vtx  = [nds{i}(I,:) ; ctr{ip2}(I,:) ; ctr{ip1}(I,:)];
    elt  = reshape((1:3*length(I))',length(I),3);
    tmp  = msh(vtx,elt,col(I));
    mesh = union(mesh,tmp);
    
    % Second elements 
    vtx  = [nds{ip1}(I,:) ; nds{ip2}(I,:) ; ctr{ip1}(I,:)];
    elt  = reshape((1:3*length(I))',length(I),3);
    tmp  = msh(vtx,elt,col(I));
    mesh = union(mesh,tmp);
    
    % Third elements
    vtx  = [nds{ip1}(I,:) ; ctr{ip1}(I,:) ; ctr{ip2}(I,:) ; ];
    elt  = reshape((1:3*length(I))',length(I),3);
    tmp  = msh(vtx,elt,col(I));
    mesh = union(mesh,tmp);
end

% Security
if size(mesh.elt,1) ~= 3*Nelt
    error('mshMidpoint2.m : unavailable case 1')
end
if size(mesh.vtx,1) ~= Nvtx+2*Nelt
    error('mshMidpoint2.m : unavailable case 2')
end
end


function mesh = mshMidpoint3(mesh)    
% Mesh nodes and edges center
[nds,ctr] = data(mesh);
   
% Refined mesh initialization with centered triangle
Nelt = length(mesh);
vtx  = [ctr{1} ; ctr{2} ; ctr{3}];
elt  = reshape((1:3*Nelt)',Nelt,3);
col  = mesh.col;
mesh = msh(vtx,elt,col);

% For each node
for i = 1:3
    % Neighbours
    ip1 = mod(i,3)+1;
    ip2 = mod(ip1,3)+1;

    % New submesh with nodes triangles
    vtx = [nds{i} ; ctr{ip2} ; ctr{ip1}];
    tmp = msh(vtx,elt,col);
    
    % Mesh fusion with previous submeshes
    mesh = union(mesh,tmp);
end

% Security
if size(mesh.elt,1) ~= 4*Nelt
    error('mshMidpoint3.m : unavailable case 1')
end
end


function [nds,ctr] = data(mesh)
% Triangle nodes
nds{1} = mesh.vtx(mesh.elt(:,1),:);
nds{2} = mesh.vtx(mesh.elt(:,2),:);
nds{3} = mesh.vtx(mesh.elt(:,3),:);

% Edges center
ctr{1} = 0.5 * (nds{2} + nds{3});
ctr{2} = 0.5 * (nds{3} + nds{1});
ctr{3} = 0.5 * (nds{1} + nds{2});
end

