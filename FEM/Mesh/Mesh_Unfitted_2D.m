classdef Mesh_Unfitted_2D < Mesh_Unfitted
    
    methods
        function obj = Mesh_Unfitted_2D(fitted_mesh,fitted_geom_interpolation)
            obj@Mesh_Unfitted(fitted_mesh,fitted_geom_interpolation);
            obj.max_subcells = 6;
            obj.nnodes_subcell = 3;
        end
        
        function [P,active_nodes] = findCutPoints_Iso(obj)
            pos_nodes = obj.fitted_geom_interpolation.pos_nodes;
            
            gamma_1 = permute(obj.x_fitted(obj.fitted_mesh.connec(obj.cut_cells,:)),[2 3 1]);
            gamma_2 = permute([obj.x_fitted(obj.fitted_mesh.connec(obj.cut_cells,2:end)),obj.x_fitted(obj.fitted_mesh.connec(obj.cut_cells,1))],[2 3 1]);
            
            P1 = repmat(pos_nodes,[1 1 size(obj.cut_cells)]);
            P2 = repmat(circshift(pos_nodes,[size(pos_nodes,1)-1 0]),[1 1 size(obj.cut_cells)]);
            P = P1+gamma_1.*(P2-P1)./(gamma_1-gamma_2);
            
            active_nodes = sign(gamma_1.*gamma_2)<=0;
        end
        
        function [P,active_nodes] = findCutPoints_Global(obj)
            index1 = permute(obj.fitted_mesh.connec(obj.cut_cells,:),[2 3 1]);
            % !! DO THIS PATCH WITH circshift LIKE IN findCutPoints_Iso !!
            index2 = [permute(obj.fitted_mesh.connec(obj.cut_cells,2:end),[2 3 1]);...
                permute(obj.fitted_mesh.connec(obj.cut_cells,1),[2 3 1])];
            gamma_1 = obj.x_fitted(index1);
            gamma_2 = obj.x_fitted(index2);
            coord1 = obj.fitted_mesh.coord(:,1); coord2 = obj.fitted_mesh.coord(:,2);
            
            P1 = [coord1(index1) coord2(index1)];
            P2 = [coord1(index2) coord2(index2)];
            P = P1+gamma_1.*(P2-P1)./(gamma_1-gamma_2);
            
            active_nodes = sign(gamma_1.*gamma_2)<=0;
        end
        
        function computeDvoluCut(obj)
            x1 = obj.unfitted_cut_coord_iso_per_cell(:,1,1);   y1 = obj.unfitted_cut_coord_iso_per_cell(:,1,2);
            x2 = obj.unfitted_cut_coord_iso_per_cell(:,2,1);   y2 = obj.unfitted_cut_coord_iso_per_cell(:,2,2);
            x3 = obj.unfitted_cut_coord_iso_per_cell(:,3,1);   y3 = obj.unfitted_cut_coord_iso_per_cell(:,3,2);
            obj.dvolu_cut = 0.5*abs((x2-x1).*(y3-y1)-(x3-x1).*(y2-y1));
        end
        
        function assignUnfittedCutCoordIsoPerCell(obj,new_unfitted_cut_coord_iso,new_subcell_cut_interior_connec_iso,c0,c1)
            new_unfitted_cut_coord_iso_x = new_unfitted_cut_coord_iso(:,1);
            new_unfitted_cut_coord_iso_y = new_unfitted_cut_coord_iso(:,2);
            
            obj.unfitted_cut_coord_iso_per_cell(c0+1:c1,:,1) = new_unfitted_cut_coord_iso_x(new_subcell_cut_interior_connec_iso);
            obj.unfitted_cut_coord_iso_per_cell(c0+1:c1,:,2) = new_unfitted_cut_coord_iso_y(new_subcell_cut_interior_connec_iso);
        end
    end
end

