classdef InclusionTraining < handle

    properties (Access = public)

    end

    properties (Access = private)

    end

    properties (Access = private)
        r
        xmin
        xmax
        ymin
        ymax
        cy
        cx
        Nr
        Ntheta
        connec
        CoarseOrder
    end

    methods (Access = public)

        function obj = InclusionTraining()
            close all
            obj.init()

            for i=1:length(obj.r)
                m = obj.mesh_rectangle_via_triangles(obj.r(i));
%                  m = obj.meshFromPython(obj.r(i));
                if i==1
                    obj.connec = m.connec;
                end

                if norm(m.connec-obj.connec)>0
                    break
                end

                %                 figure
                %                 m.plot()
                %                 m.plotAllNodes();
                %                 m = obj.createReferenceMesh();
                data = Training(m,obj.CoarseOrder);
%                                 obj.printdisplacements(data.uSbd,m,i)
                p = OfflineDataProcessor(data);
                EIFEoper = p.computeROMbasis();
                EIFEoper.U = EIFEoper.Udef + EIFEoper.Urb;
                EIFEoper.U = EIFEoper.U(:);
%                 EIFEoper.Kfine = data.LHSsbd;
                EIFEoper.snapshots = data.uSbd;
%                 filePath = ['/home/raul/Documents/GitHub/EPFL/test/data_' num2str(obj.r(i), '%.3f') '.mat'];
                filePath = ['./EPFL/dataSquare/data_' num2str(obj.r(i), '%.3f') '.mat'];
                save(filePath,'EIFEoper')
            end
        end

    end

    methods (Access = private)

        function init(obj)
%             N = 80;
%             % Interval bounds
%             a = 0.99;
%             b = 0.8;
%             % Index vector
%             i = 0:N;
%             % Cosine spacing formula
%             obj.r = (a + b)/2 + (b - a)/2 * cos(pi * (1 - i / N));
            obj.r    = 0.8:0.01:0.961;
            obj.xmin = -1;
            obj.xmax = 1;
            obj.ymin = -1;
            obj.ymax = 1;
            obj.cx = 0;
            obj.cy = 0;
            obj.Nr=7;
            obj.Ntheta=14;
            obj.CoarseOrder = 1;
        end

        function printdisplacements(obj,Usbd,mesh,ind)
            for i = 1:size(Usbd,2)
                v = Usbd(:,i);
                EIFEMtesting.plotSolution(v,mesh,1,ind,i,[])
            end
        end

        function mesh = createReferenceMesh(obj)

            %UnitMesh better
            %             x1      = linspace(-1,1,50);
            %             x2      = linspace(-1,1,50);
            %             [xv,yv] = meshgrid(x1,x2);
            %             [F,V]   = mesh2tri(xv,yv,zeros(size(xv)),'x');
            %             s.coord  = V(:,1:2);
            %              s.coord(s.coord(:,1)== obj.xmax & s.coord(:,2)==obj.ymax,:) =...
            %                 s.coord(s.coord(:,1)== obj.xmax & s.coord(:,2)==obj.ymax,:)-[1e-9,0];
            %             s.coord(s.coord(:,1)== obj.xmax & s.coord(:,2)==obj.ymin,:) =...
            %                 s.coord(s.coord(:,1)== obj.xmax & s.coord(:,2)==obj.ymin,:)-[1e-9,0];
            %             s.coord(s.coord(:,1)== obj.xmin & s.coord(:,2)==obj.ymax,:) =...
            %                 s.coord(s.coord(:,1)== obj.xmin & s.coord(:,2)==obj.ymax,:)+[1e-9,0];
            %             s.coord(s.coord(:,1)== obj.xmin & s.coord(:,2)==obj.ymin,:) =...
            %                 s.coord(s.coord(:,1)== obj.xmin & s.coord(:,2)==obj.ymin,:)+[1e-9,0];
            %             s.connec = F;
            %             mesh = Mesh.create(s);

            mesh = QuadMesh(2,2,50,50);
            delta= 1e-9;
            s.coord = mesh.coord-[1,1];
            s.coord(s.coord(:,1)== obj.xmax & s.coord(:,2)==obj.ymax,:) =...
                s.coord(s.coord(:,1)== obj.xmax & s.coord(:,2)==obj.ymax,:)+[-delta,-delta];
            s.coord(s.coord(:,1)== obj.xmax & s.coord(:,2)==obj.ymin,:) =...
                s.coord(s.coord(:,1)== obj.xmax & s.coord(:,2)==obj.ymin,:)+[-delta,+delta];
            s.coord(s.coord(:,1)== obj.xmin & s.coord(:,2)==obj.ymax,:) =...
                s.coord(s.coord(:,1)== obj.xmin & s.coord(:,2)==obj.ymax,:)+[+delta,-delta];
            s.coord(s.coord(:,1)== obj.xmin & s.coord(:,2)==obj.ymin,:) =...
                s.coord(s.coord(:,1)== obj.xmin & s.coord(:,2)==obj.ymin,:)+[+delta,+delta];
            s.connec=mesh.connec;
            mesh = Mesh.create(s);

        end

        function [nodes, elements] = mesh_triangle_sector(obj,cx, cy, corner1, corner2, r_inner, Nr, Ntheta)
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
%                 r_inner = obj.ray_square_intersection(dir, r_inner);
                %uncomment for square
                 r_inner_dir = r_inner;
%                 r_inner_dir = r_inner / max(abs(dir)); 

                % Create points from r_inner to r_max along dir
                r_vals = linspace(r_inner_dir, r_max, Nr + 1)';

                %for circular inclusion
%                 r_vals = linspace(r_inner, r_max, Nr + 1)';
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

        function r_inner = ray_square_intersection(obj,dir, a)
    dx = dir(1);
    dy = dir(2);

    % Very small tolerance to avoid division issues
    epsv = 1e-12;

    % Compute intersection distances to all 4 sides
    t = [];

    if abs(dx) > epsv
        t1 = ( a) / dx;
        t2 = (-a) / dx;
        t = [t, t1, t2];
    end

    if abs(dy) > epsv
        t3 = ( a) / dy;
        t4 = (-a) / dy;
        t = [t, t3, t4];
    end

    % Keep only intersections in front of the ray (t > 0)
    t = t(t > 0);

    % Compute the hit coordinates for all candidate intersections
    hits = [t' * dx, t' * dy];

    % Keep only hits lying on the square boundary
    valid = abs(hits(:,1)) <= a + 1e-10 & abs(hits(:,2)) <= a + 1e-10;

    t = t(valid);

    % Choose the smallest positive valid intersection
    r_inner = min(t);
        end


        function mesh = mesh_rectangle_via_triangles(obj,r)
            % Mesh rectangle with hole by meshing 4 triangles and flipping the mesh

            % Define corners
            corners = [obj.xmax obj.ymax;
                obj.xmin obj.ymax;
                obj.xmin obj.ymin;
                obj.xmax obj.ymin];

            nodes_all = [];
            elements_all = [];
            node_offset = 0;

            % Loop over 4 corners - mesh each triangle sector and transform it
            for k = 1:4
                c1 = corners(k,:);
                c2 = corners(mod(k,4)+1,:);

                % Mesh one triangle sector
                [nodes, elements] = obj.mesh_triangle_sector(obj.cx, obj.cy, c1, c2, r, obj.Nr, obj.Ntheta);

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
            delta= 1e-9;
            s.coord(s.coord(:,1)== obj.xmax & s.coord(:,2)==obj.ymax,:) =...
                s.coord(s.coord(:,1)== obj.xmax & s.coord(:,2)==obj.ymax,:)+[-delta,-delta];
            s.coord(s.coord(:,1)== obj.xmax & s.coord(:,2)==obj.ymin,:) =...
                s.coord(s.coord(:,1)== obj.xmax & s.coord(:,2)==obj.ymin,:)+[-delta,+delta];
            s.coord(s.coord(:,1)== obj.xmin & s.coord(:,2)==obj.ymax,:) =...
                s.coord(s.coord(:,1)== obj.xmin & s.coord(:,2)==obj.ymax,:)+[+delta,-delta];
            s.coord(s.coord(:,1)== obj.xmin & s.coord(:,2)==obj.ymin,:) =...
                s.coord(s.coord(:,1)== obj.xmin & s.coord(:,2)==obj.ymin,:)+[+delta,+delta];
%             s.coord(s.coord(:,1)== obj.xmax & s.coord(:,2)==obj.ymax,:) =...
%                 s.coord(s.coord(:,1)== obj.xmax & s.coord(:,2)==obj.ymax,:)-[1e-9,0];
%             s.coord(s.coord(:,1)== obj.xmax & s.coord(:,2)==obj.ymin,:) =...
%                 s.coord(s.coord(:,1)== obj.xmax & s.coord(:,2)==obj.ymin,:)-[1e-9,0];
%             s.coord(s.coord(:,1)== obj.xmin & s.coord(:,2)==obj.ymax,:) =...
%                 s.coord(s.coord(:,1)== obj.xmin & s.coord(:,2)==obj.ymax,:)+[1e-9,0];
%             s.coord(s.coord(:,1)== obj.xmin & s.coord(:,2)==obj.ymin,:) =...
%                 s.coord(s.coord(:,1)== obj.xmin & s.coord(:,2)==obj.ymin,:)+[1e-9,0];

            %             if isempty(obj.connec)
            %                 s.connec = elements_final;
            %                 obj.connec = s.connec;
            %             else
            %                 s.connec = obj.connec;
            %             end

            s.connec = elements_final;
            mesh = Mesh.create(s);
        end

        function m = meshFromPython(obj,r)
            mesh_module= py.importlib.import_module('mesh_module');

            % Inputs
            n = py.list({int32(obj.Nr), int32(obj.Ntheta)});
            degree = int32(1);
            p0 = py.list({obj.xmin, obj.ymin});
            p1 = py.list({obj.xmax, obj.ymax});
            parameters = py.numpy.array([r, r, r, r]);

            % Call Python function
            res = mesh_module.create_subdomain_mesh(n, degree, p0, p1, parameters);
            coords = res{1};
            s.coord = double(coords);
            tol = 1e-5;

            mask = abs(s.coord(:,1) - obj.xmax) < tol;
            s.coord(mask, 1) = obj.xmax;

             mask = abs(s.coord(:,1) - obj.xmin) < tol;
            s.coord(mask, 1) = obj.xmin;

             mask = abs(s.coord(:,2) - obj.ymax) < tol;
            s.coord(mask, 2) = obj.ymax;

             mask = abs(s.coord(:,2) - obj.ymin) < tol;
            s.coord(mask, 2) = obj.ymin;

            % Top-right corner (xmax, ymax)
            mask = abs(s.coord(:,1) - obj.xmax) < tol & abs(s.coord(:,2) - obj.ymax) < tol;
            s.coord(mask, :) = s.coord(mask, :) - [1e-9, 0];

            % Bottom-right corner (xmax, ymin)
            mask = abs(s.coord(:,1) - obj.xmax) < tol & abs(s.coord(:,2) - obj.ymin) < tol;
            s.coord(mask, :) = s.coord(mask, :) - [1e-9, 0];

            % Top-left corner (xmin, ymax)
            mask = abs(s.coord(:,1) - obj.xmin) < tol & abs(s.coord(:,2) - obj.ymax) < tol;
            s.coord(mask, :) = s.coord(mask, :) + [1e-9, 0];

            % Bottom-left corner (xmin, ymin)
            mask = abs(s.coord(:,1) - obj.xmin) < tol & abs(s.coord(:,2) - obj.ymin) < tol;
            s.coord(mask, :) = s.coord(mask, :) + [1e-9, 0];

            conn = res{2};
            conn_cell = cell(conn);
            conn_mat = cellfun(@(x) double(x), conn_cell, 'UniformOutput', false);
            conn_array = cell2mat(conn_mat');  % Transpose first to get elements in rows
            conn_array(:, [3, end]) = conn_array(:, [end, 3]);
            conn_array = conn_array + 1;

            s.connec = conn_array;
            m = Mesh.create(s);
        end
    end

end