classdef Mesh_Unfitted_2D < Mesh_Unfitted
    properties (GetAccess = public, SetAccess = protected)
        ndim = 2;
    end
    
    methods (Access = public)
        function [P,active_nodes] = computeCutPoints_Iso(obj)
            pos_nodes = obj.fitted_geom_interpolation.pos_nodes;
            
            gamma_1 = permute(obj.x_fitted(obj.fitted_mesh.connec(obj.fitted_cut_cells,:)),[2 3 1]);
            gamma_2 = circshift(gamma_1,[-1 0 0]);
            
            P1 = repmat(pos_nodes,[1 1 size(obj.fitted_cut_cells)]);
            P2 = circshift(P1,[-1 0 0]);
            P = P1 + gamma_1.*(P2-P1)./(gamma_1-gamma_2);
            
            active_nodes = sign(gamma_1.*gamma_2)<=0;
        end
        
        function [P,active_nodes] = computeCutPoints_Global(obj)
            index1 = permute(obj.fitted_mesh.connec(obj.fitted_cut_cells,:),[2 3 1]);
            index2 = circshift(index1,[-1 0 0]);
            
            gamma_1 = obj.x_fitted(index1);
            gamma_2 = obj.x_fitted(index2);
           
            coord1 = obj.fitted_mesh.coord(:,1);
            coord2 = obj.fitted_mesh.coord(:,2);
            
            P1 = [coord1(index1) coord2(index1)];
            P2 = [coord1(index2) coord2(index2)];
            P = P1 + gamma_1.*(P2-P1)./(gamma_1-gamma_2);
            
            active_nodes = sign(gamma_1.*gamma_2)<=0;
        end
        
        function assignUnfittedCutCoordIsoPerCell(obj,new_coord_iso,new_interior_subcell_connec,lowerBound,upperBound)
            new_coord_iso_x = new_coord_iso(:,1);
            new_coord_iso_y = new_coord_iso(:,2);
            
            obj.coord_iso_per_cell(lowerBound+1:upperBound,:,1) = new_coord_iso_x(new_interior_subcell_connec);
            obj.coord_iso_per_cell(lowerBound+1:upperBound,:,2) = new_coord_iso_y(new_interior_subcell_connec);
        end
    end
end

