classdef Mesh_Unfitted_3D < Mesh_Unfitted
    methods (Access = public)
        function obj = Mesh_Unfitted_3D
            obj.ndim = 3;
        end

%     end
%     methods (Access = private)
        function [P,active_nodes] = computeCutPoints_Iso(obj)
            pos_nodes = obj.background_geom_interpolation.pos_nodes;
            
            iteration_1 = obj.background_geom_interpolation.iteration(1,:);
            iteration_2 = obj.background_geom_interpolation.iteration(2,:);
            
            gamma_1 = permute(obj.x_background(obj.background_mesh.connec(obj.background_cut_cells,iteration_1)),[2 3 1]);
            gamma_2 = permute(obj.x_background(obj.background_mesh.connec(obj.background_cut_cells,iteration_2)),[2 3 1]);
            
            P1 = repmat(pos_nodes(iteration_1,:),[1 1 size(obj.background_cut_cells)]);
            P2 = repmat(pos_nodes(iteration_2,:),[1 1 size(obj.background_cut_cells)]);
            P = P1+gamma_1.*(P2-P1)./(gamma_1-gamma_2);
            
            active_nodes = sign(gamma_1.*gamma_2)<=0;
        end
        
        function [P,active_nodes] = computeCutPoints_Global(obj)
            iteration_1 = obj.background_geom_interpolation.iteration(1,:);
            iteration_2 = obj.background_geom_interpolation.iteration(2,:);
            
            index1 = permute(obj.background_mesh.connec(obj.background_cut_cells,iteration_1),[2 3 1]);
            index2 = permute(obj.background_mesh.connec(obj.background_cut_cells,iteration_2),[2 3 1]);
           
            gamma_1 = obj.x_background(index1);
            gamma_2 = obj.x_background(index2);
            
            coord1 = obj.background_mesh.coord(:,1);
            coord2 = obj.background_mesh.coord(:,2);
            coord3 = obj.background_mesh.coord(:,3);
            
            P1 = [coord1(index1) coord2(index1) coord3(index1)];
            P2 = [coord1(index2) coord2(index2) coord3(index2)];
            P = P1+gamma_1.*(P2-P1)./(gamma_1-gamma_2);
            
            active_nodes = sign(gamma_1.*gamma_2)<=0;
        end
        
        function assignUnfittedCutCoordIsoPerCell(obj,new_coord_iso,new_interior_subcell_connec,lowerBound,upperBound)
            new_coord_iso_x = new_coord_iso(:,1);
            new_coord_iso_y = new_coord_iso(:,2);
            new_coord_iso_z = new_coord_iso(:,3);
            
            obj.coord_iso_per_cell(lowerBound+1:upperBound,:,1) = new_coord_iso_x(new_interior_subcell_connec);
            obj.coord_iso_per_cell(lowerBound+1:upperBound,:,2) = new_coord_iso_y(new_interior_subcell_connec);
            obj.coord_iso_per_cell(lowerBound+1:upperBound,:,3) = new_coord_iso_z(new_interior_subcell_connec);
        end
    end
end

