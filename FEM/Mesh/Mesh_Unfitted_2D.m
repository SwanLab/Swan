classdef Mesh_Unfitted_2D < Mesh_Unfitted
    
    methods
        function obj = Mesh_Unfitted_2D(fitted_mesh,x_fitted,fitted_geom_interpolation)
            obj@Mesh_Unfitted(fitted_mesh,fitted_geom_interpolation,x_fitted);
            obj.max_subcells = 6;
            obj.nnodes_subcell = 3;
        end
        
        function [subcell_cut_interior_coord_iso,subcell_cut_interior_coord_global,subcell_cut_interior_x_value,subcell_cut_interior_connec_iso] = computeSubcells(obj,fitted_connec,subcell_cutPoints_iso,subcell_cutPoints_global)
            [subcell_cut_interior_coord_iso,subcell_cut_interior_coord_global,subcell_cut_interior_x_value,subcell_cut_interior_connec_iso] = obj.computeInteriorSubcells(fitted_connec,subcell_cutPoints_iso,subcell_cutPoints_global);
        end
        
        function [subcell_cut_interior_coord_iso,subcell_cut_interior_coord_global,subcell_cut_interior_x_value,subcell_cut_interior_connec] = computeInteriorSubcells(obj,fitted_cell_connec,subcell_cutPoints_iso,subcell_cutPoints_global)
            subcell_coord_iso = [obj.fitted_geom_interpolation.pos_nodes; subcell_cutPoints_iso];
            subcell_coord_global = [obj.fitted_mesh.coord(fitted_cell_connec,:); subcell_cutPoints_global];
            subcell_x_value = [obj.x_fitted(fitted_cell_connec); zeros(size(subcell_cutPoints_iso,1),1)]';
            
            subcell_cut_interior_connec = obj.computeInteriorSubcellsConnectivities(subcell_coord_iso,subcell_x_value);
            
%             interior_nodes = unique(subcell_cut_interior_connec);
%             subcell_cut_interior_coord_iso = subcell_coord_iso(interior_nodes,:);
%             subcell_cut_interior_coord_global = subcell_coord_global(interior_nodes,:);
%             subcell_cut_interior_x_value = subcell_x_value(interior_nodes);
            
            subcell_cut_interior_coord_iso = subcell_coord_iso;
            subcell_cut_interior_coord_global = subcell_coord_global;
            subcell_cut_interior_x_value = subcell_x_value;
end
        
        function [P,active_nodes] = findCutPoints_Iso(obj)
            pos_nodes = obj.fitted_geom_interpolation.pos_nodes;
            
            gamma_1 = permute(obj.x_fitted(obj.fitted_mesh.connec(obj.cut_cells,:)),[2 3 1]);
            gamma_2 = permute([obj.x_fitted(obj.fitted_mesh.connec(obj.cut_cells,2:end)),obj.x_fitted(obj.fitted_mesh.connec(obj.cut_cells,1))],[2 3 1]);
            
            P1 = repmat(pos_nodes,[1 1 size(obj.cut_cells)]);
            % !! There's an specific Matlab function to do this !!
            P2 = repmat([pos_nodes(2:end,:);pos_nodes(1,:)],[1 1 size(obj.cut_cells)]);
            P = P1+gamma_1.*(P2-P1)./(gamma_1-gamma_2);
            
            active_nodes = sign(gamma_1.*gamma_2)<0;
        end
        
        function [P,active_nodes] = findCutPoints_Global(obj)
            index1 = permute(obj.fitted_mesh.connec(obj.cut_cells,:),[2 3 1]);
            index2 = [permute(obj.fitted_mesh.connec(obj.cut_cells,2:end),[2 3 1]);...
                permute(obj.fitted_mesh.connec(obj.cut_cells,1),[2 3 1])];
            gamma_1 = obj.x_fitted(index1);
            gamma_2 = obj.x_fitted(index2);
            coord1 = obj.fitted_mesh.coord(:,1); coord2 = obj.fitted_mesh.coord(:,2);
            
            P1 = [coord1(index1) coord2(index1)];
            P2 = [coord1(index2) coord2(index2)];
            P = P1+gamma_1.*(P2-P1)./(gamma_1-gamma_2);
            
            active_nodes = sign(gamma_1.*gamma_2)<0;
        end
        
        function computeDvoluCut(obj)
%             x1 = obj.unfitted_cut_coord_iso(obj.unfitted_cut_connec_iso(:,1),1); y1 = obj.unfitted_cut_coord_iso(obj.unfitted_cut_connec_iso(:,1),2);
%             x2 = obj.unfitted_cut_coord_iso(obj.unfitted_cut_connec_iso(:,2),1); y2 = obj.unfitted_cut_coord_iso(obj.unfitted_cut_connec_iso(:,2),2);
%             x3 = obj.unfitted_cut_coord_iso(obj.unfitted_cut_connec_iso(:,3),1); y3 = obj.unfitted_cut_coord_iso(obj.unfitted_cut_connec_iso(:,3),2);
%             
%             obj.dvolu_cut= 0.5*abs((x2-x1).*(y3-y1)-(x3-x1).*(y2-y1));
            
            x1 = obj.unfitted_cut_coord_iso_per_cell(:,1,1); y1 = obj.unfitted_cut_coord_iso_per_cell(:,1,2); x2 = obj.unfitted_cut_coord_iso_per_cell(:,2,1);
            y2 = obj.unfitted_cut_coord_iso_per_cell(:,2,2); x3 = obj.unfitted_cut_coord_iso_per_cell(:,3,1); y3 = obj.unfitted_cut_coord_iso_per_cell(:,3,2);
            obj.dvolu_cut = 0.5*abs((x2-x1).*(y3-y1)-(x3-x1).*(y2-y1));
        end
    end
    
    methods (Static)
        function subcell_connec = computeInteriorSubcellsConnectivities(subcell_coord_iso,subcell_x_value)
            DT = delaunayTriangulation(subcell_coord_iso);
            subcell_connec = DT.ConnectivityList;
            is_interior = all(subcell_x_value(subcell_connec) <= 0,2);
            
            subcell_connec = subcell_connec(is_interior,:);
        end
        
        function connec_global = computeGlobalConnectivities(coord_global)
            DT = delaunayTriangulation(coord_global);
            connec_global = DT.ConnectivityList;
        end
    end
end

