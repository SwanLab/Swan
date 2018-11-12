function [corners,pnods] = compute_periodic_boundary_nodes_square(coordinates)
% PERIODIC BOUNDARY COND
% creation of a list containing the couples that define the periodicity
% 1) list of nodes on each side (two vertical, two horizontal)
% 2) sort each list based on the corresponding coordinate
% 3) generation of the couples and store them in 'pnods'
% 
% remark: 
% [Y,I] = SORT(X,DIM,MODE) also returns an index matrix I.
% If X is a vector, then Y = X(I).  

% nodes in the left-vertical side, without the corners
href = 0.025;
h=href; L=[];Y=[];
for i=1:size(coordinates,1)
    if (coordinates(i,1)<h/5 && coordinates(i,2)>h/5 && coordinates(i,2)<1-h/5 )
        L = [L i];
        Y = [Y coordinates(i,2)];
    end
end
[Y1,I] = sort(Y);
V1 = L(I);

% nodes in the right-vertical side, without the corners
h=href; L=[];Y=[];
for i=1:size(coordinates,1)
    if (coordinates(i,1)>1-h/5 && coordinates(i,2)>h/5 && coordinates(i,2)<1-h/5)
        L = [L i];
        Y = [Y coordinates(i,2)];
    end
end
[Y1,I] = sort(Y);
V2 = L(I);

% nodes in the bottom-horizontal side, without the corners
h=href; L=[];X=[];
for i=1:size(coordinates,1)
    if (coordinates(i,2)<h/5 && coordinates(i,1)>h/5 && coordinates(i,1)<1-h/5 )
        L = [L i];
        X = [X coordinates(i,1)];
    end
end
[X1,I] = sort(X);
H1 = L(I);

% nodes in the top-horizontal side, without the corners
h=href; L=[];X=[];
for i=1:size(coordinates,1)
    if (coordinates(i,2)>1-h/5 && coordinates(i,1)>h/5 && coordinates(i,1)<1-h/5)
        L = [L i];
        X = [X coordinates(i,1)];
    end
end
[X1,I] = sort(X);
H2 = L(I);
pnods = [V1 H1; V2 H2]; % lista de nodos 

nodes = 1:length(coordinates(:,1));
%Vertex 1 y 2
index_bin_max_y = coordinates(:,2) == max(coordinates(:,2));
index_max_y = nodes(index_bin_max_y);

[~,index1_local] = min(coordinates(index_max_y,1));
index_vertex_1 = index_max_y(index1_local);
[~,index2_local] = max(coordinates(index_max_y,1));
index_vertex_2 = index_max_y(index2_local);

%Vertex 3 y 4
index_bin_min_y = coordinates(:,2) == min(coordinates(:,2));
index_min_y = nodes(index_bin_min_y);

[~,index3_local] = max(coordinates(index_min_y,1));
index_vertex_3 = index_min_y(index3_local);
[~,index4_local] = min(coordinates(index_min_y,1));
index_vertex_4 = index_min_y(index4_local);

index_vertex = [index_vertex_1;index_vertex_2;index_vertex_3;index_vertex_4];
corners = index_vertex';


end
% element.lglib = [2*1-1    2*3277-1  2*1     2*3277 ;  ...
%                  2*4225-1 2*3278-1  2*4225  2*3278];
% element.l2penalty = [   1  3277    1 3278;
%                      3278  4225 3277 4225];
% element.penalty = 10000000*1*0.015625^2;

%  1    3277
%
% 3278  4225

