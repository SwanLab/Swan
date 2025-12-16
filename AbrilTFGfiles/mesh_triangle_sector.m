function [nodes, elements] = mesh_triangle_sector(cx, cy, corner1, corner2, r_inner, Nr, Ntheta)
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
