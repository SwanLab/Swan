function Ne=vertex_neighbours(FV)
% This function VERTEX_NEIGHBOURS will search in a face list for all 
% the neigbours of each vertex.
%
% Ne=vertex_neighbours(FV)
%

Ne=vertex_neighbours_double(FV.faces(:,1),FV.faces(:,2),FV.faces(:,3),FV.vertices(:,1),FV.vertices(:,2),FV.vertices(:,3));
