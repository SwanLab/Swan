function mesh = mesh_rectangle_via_triangles(r,xmax,xmin,ymax,ymin,Nr,Ntheta,x0,y0)
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
        [nodes, elements] = mesh_triangle_sector(x0, y0, c1, c2, r, Nr, Ntheta);
        % Append nodes, elements with offset
        elements = elements+node_offset;
        nodes_all = [nodes_all; nodes];
        elements_all = [elements_all; elements];
        node_offset = size(nodes_all,1);
    end
    tol = 1e-8;
    % Step 1: Snap coordinates to a regular grid
    rounded_nodes = round(nodes_all / tol) * tol;
    %             rounded_nodes = round(nodes_all * 1e12) / 1e12;
    [~, ia, ic] = unique(rounded_nodes, 'rows', 'stable');
    nodes_final = nodes_all(ia, :);
    elements_final = ic(elements_all);
    %             [~, ~, J] = uniquetol(nodes_all, 1e-12, 'ByRows', true);
    %             firstidx = accumarray(J, (1:length(J))', [], @min);
    %             nodes_final = nodes_all(firstidx, :);
    %             elements_final = J(elements_all);
    %
    % %             % Merge duplicate nodes (within a tolerance)
    % %             [nodes_final, ~, ic] = uniquetol(nodes_all, 1e-12, 'ByRows', true,'Stable');
    % %
    % %             % Re-map element indices to unique node list
    % %             elements_final = ic(elements_all);
    s.coord = nodes_final;
    s.coord(s.coord(:,1)== xmax & s.coord(:,2)==ymax,:) =...
        s.coord(s.coord(:,1)== xmax & s.coord(:,2)==ymax,:)+[-1e-9,-1e-9];
    s.coord(s.coord(:,1)== xmax & s.coord(:,2)==ymin,:) =...
        s.coord(s.coord(:,1)== xmax & s.coord(:,2)==ymin,:)+[-1e-9,1e-9];
    s.coord(s.coord(:,1)== xmin & s.coord(:,2)==ymax,:) =...
        s.coord(s.coord(:,1)== xmin & s.coord(:,2)==ymax,:)+[1e-9,-1e-9];
    s.coord(s.coord(:,1)== xmin & s.coord(:,2)==ymin,:) =...
        s.coord(s.coord(:,1)== xmin & s.coord(:,2)==ymin,:)+[1e-9,1e-9];
    %             if isempty(obj.connec)
    %                 s.connec = elements_final;
    %                 obj.connec = s.connec;
    %             else
    %                 s.connec = obj.connec;
    %             end
    s.connec = elements_final;
    mesh = Mesh.create(s);
end        