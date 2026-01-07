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
                elements = elements + node_offset;
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
%                 r_inner_dir = r_inner; % for circular inclusion
                 r_inner_dir = r_inner / max(abs(dir)); % for sqaure
                % Create points from r_inner to r_max along dir
                r_vals = linspace(r_inner_dir, r_max, Nr + 1)';
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
