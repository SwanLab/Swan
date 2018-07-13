function [corners,pnods] = compute_periodic_boundary_nodes_hexagonal(coordinates)
% PERIODIC BOUNDARY COND
% creation of a list containing the couples that define the periodicity
% 1) list of nodes on each side (two vertical, two horizontal)
% 2) sort each list based on the corresponding coordinate
% 3) generation of the couples and store them in 'pnods'
href = 0.025; 

nodes = 1:length(coordinates(:,1));
%Vertex 1 
[~,index_vertex_1] = min(coordinates(:,1));
%Vertex 2 y 3
index_bin_max_y = coordinates(:,2) == max(coordinates(:,2));
index_max_y = nodes(index_bin_max_y);

[~,index2_local] = min(coordinates(index_max_y,1));
index_vertex_2 = index_max_y(index2_local);
[~,index3_local] = max(coordinates(index_max_y,1));
index_vertex_3 = index_max_y(index3_local);

%Vertex  4
[~,index_vertex_4] = max(coordinates(:,1));

%Vertex 5 y 6
index_bin_min_y = coordinates(:,2) == min(coordinates(:,2));
index_min_y = nodes(index_bin_min_y);

[~,index5_local] = max(coordinates(index_min_y,1));
index_vertex_5 = index_min_y(index5_local);
[~,index6_local] = min(coordinates(index_min_y,1));
index_vertex_6 = index_min_y(index6_local);

index_vertex = [index_vertex_1;index_vertex_2;index_vertex_3;index_vertex_4;index_vertex_5;index_vertex_6];

nsides = 6;
order_corners = [1:nsides,1];
for isides = 1:nsides
    idx1 = order_corners(isides);
    idx2 = order_corners(isides+1);
    
    line_vector = coordinates(index_vertex(idx2),:) - coordinates(index_vertex(idx1),:);
    point_line = coordinates(index_vertex(idx1),:);
    
    distance_point2line = @(x,y) (line_vector(1)*(y-point_line(2)) - line_vector(2)*(x-point_line(1)))/norm(line_vector);
    dist2line = distance_point2line(coordinates(:,1),coordinates(:,2));
    
    index_nodes(:,isides) = setdiff(nodes(abs(dist2line) < href),index_vertex); % nodes in the reference line
end


master = index_nodes(:,1:3);
slave = index_nodes(:,4:6);

pnods = [master(:) slave(:)]';
corners = index_vertex';


end



