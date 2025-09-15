xmin = -1; xmax = 1;
ymin = -1; ymax = 1;
cx = 0; cy = 0;
r = 0.3;

Nr = 15; % radial divisions
Ntheta = 30; % angular divisions per triangle sector

[nodes, elements] = mesh_rectangle_via_triangles(cx, cy, r, xmin, xmax, ymin, ymax, Nr, Ntheta);

figure; hold on; axis equal;
for i = 1:size(elements,1)
    elem_nodes = elements(i,:);
    patch('Vertices', nodes, 'Faces', elem_nodes, 'FaceColor', 'none', 'EdgeColor', 'k');
end

theta = linspace(0, 2*pi, 200);
plot(cx + r*cos(theta), cy + r*sin(theta), 'r-', 'LineWidth', 2);

plot([xmin xmax xmax xmin xmin], [ymin ymin ymax ymax ymin], 'b-', 'LineWidth', 2);

title('Radial mesh via 4 triangle sectors');


function [nodes, elements] = mesh_triangle_sector(cx, cy, corner1, corner2, r_inner, Nr, Ntheta)
    % Mesh radial elements in triangle formed by center and two corners
    % Inputs:
    % cx, cy - center of circle/hole
    % r - hole radius (inner radius)
    % corner1, corner2 - two rectangle corners (2x1 or 1x2 vectors)
    % Nr - number of radial divisions
    % Ntheta - number of angular divisions between corner1 and corner2
    
   v1 = corner1 - [cx, cy];
    v2 = corner2 - [cx, cy];

    % Prepare storage
    nodes = [];
    for i = 0:Ntheta
        % Linear interpolation between v1 and v2
        t = i / Ntheta;
        edge_vec = (1 - t) * v1 + t * v2;

        % Max radius along this direction
        r_max = norm(edge_vec);

        % Normalize direction vector
        dir = edge_vec / r_max;

        % Create points from r_inner to r_max along dir
        r_vals = linspace(r_inner, r_max, Nr + 1)';
        pts = [cx, cy] + r_vals * dir;

        nodes = [nodes; pts];
    end

    % Reshape node indices into a grid: (Nr+1) rows × (Ntheta+1) cols
    node_grid = reshape(1:size(nodes,1), Nr+1, Ntheta+1);

    % Build quad elements
    elements = [];
    for j = 1:Ntheta
        for i = 1:Nr
            n1 = node_grid(i, j);
            n2 = node_grid(i+1, j);
            n3 = node_grid(i+1, j+1);
            n4 = node_grid(i, j+1);
            elements(end+1, :) = [n1 n2 n3 n4];
        end
    end
end


function [nodes_all, elements_all] = mesh_rectangle_via_triangles(cx, cy, r, xmin, xmax, ymin, ymax, Nr, Ntheta)
    % Mesh rectangle with hole by meshing 4 triangles and flipping the mesh
    
    % Define corners
    corners = [xmax ymax;
               xmin ymax;
               xmin ymin;
               xmax ymin];
    
    nodes_all = [];
    elements_all = [];
    node_offset = 0;
    
    % Loop over 4 corners - mesh each triangle sector and transform it
    for k = 1:4
        c1 = corners(k,:);
        c2 = corners(mod(k,4)+1,:);
        
        % Mesh one triangle sector
        [nodes, elements] = mesh_triangle_sector(cx, cy, c1, c2,r, Nr, Ntheta);
        
        % Append nodes, elements with offset
        elements = elements + node_offset;
        nodes_all = [nodes_all; nodes];
        elements_all = [elements_all; elements];
        node_offset = size(nodes_all,1);
    end

     % Merge duplicate nodes (within a tolerance)
    [nodes_final, ~, ic] = uniquetol(nodes_all, 1e-12, 'ByRows', true);

    % Re-map element indices to unique node list
    elements_final = ic(elements_all);
end