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
    end
    
    methods (Access = public)
        
        function obj = InclusionTraining()
            obj.init()
            m = obj.mesh_rectangle_via_triangles();
            data = Training(m);
            p = OfflineDataProcessor(data);
            EIFEoper = p.computeROMbasis();
            EIFEoper.U = EIFEoper.Udef + EIFEoper.Urb;
            EIFEoper.U = EIFEoper.U(:);
            filePath = ['/home/raul/Documents/GitHub/EPFL/data_' num2str(obj.r) '.mat'];
            save(filePath,'EIFEoper')
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.r    = 0.83;
            obj.xmin = -1; 
            obj.xmax = 1;
            obj.ymin = -1;
            obj.ymax = 1;
            obj.cx = 0; 
            obj.cy = 0;
            obj.Nr=7;
            obj.Ntheta=14;
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
        
        
        function mesh = mesh_rectangle_via_triangles(obj)
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
                [nodes, elements] = obj.mesh_triangle_sector(obj.cx, obj.cy, c1, c2,obj.r, obj.Nr, obj.Ntheta);
                
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
            
            s.coord = nodes_final;
            s.coord(s.coord(:,1)== obj.xmax & s.coord(:,2)==obj.ymax,:) =...
                s.coord(s.coord(:,1)== obj.xmax & s.coord(:,2)==obj.ymax,:)-[1e-9,0];
            s.coord(s.coord(:,1)== obj.xmax & s.coord(:,2)==obj.ymin,:) =...
                s.coord(s.coord(:,1)== obj.xmax & s.coord(:,2)==obj.ymin,:)-[1e-9,0];
            s.coord(s.coord(:,1)== obj.xmin & s.coord(:,2)==obj.ymax,:) =...
                s.coord(s.coord(:,1)== obj.xmin & s.coord(:,2)==obj.ymax,:)+[1e-9,0];
            s.coord(s.coord(:,1)== obj.xmin & s.coord(:,2)==obj.ymin,:) =...
                s.coord(s.coord(:,1)== obj.xmin & s.coord(:,2)==obj.ymin,:)+[1e-9,0];


            s.connec = elements_final;
            mesh = Mesh.create(s);
        end
        
    end
    
end