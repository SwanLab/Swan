classdef Mesh_Unfitted_3D < Mesh_Unfitted
    properties (GetAccess = public, SetAccess = protected)
        ndim = 3;
    end
    
    methods (Access = public)
        function [P,active_nodes]=findCutPoints_Iso(obj)
            pos_nodes = obj.fitted_geom_interpolation.pos_nodes;
            
            iteration_1 = obj.fitted_geom_interpolation.iteration(1,:);
            iteration_2 = obj.fitted_geom_interpolation.iteration(2,:);
            
            gamma_1 = permute(obj.x_fitted(obj.fitted_mesh.connec(obj.cut_cells,iteration_1)),[2 3 1]);
            gamma_2 = permute(obj.x_fitted(obj.fitted_mesh.connec(obj.cut_cells,iteration_2)),[2 3 1]);
            
            P1 = repmat(pos_nodes(iteration_1,:),[1 1 size(obj.cut_cells)]);
            % !! There's an specific Matlab function to do this patch: circshift!!
            P2 = repmat(pos_nodes(iteration_2,:),[1 1 size(obj.cut_cells)]);
            P = P1+gamma_1.*(P2-P1)./(gamma_1-gamma_2);
            
            active_nodes = sign(gamma_1.*gamma_2)<=0;
        end
        
        function [P,active_nodes]=findCutPoints_Global(obj)
            iteration_1=obj.fitted_geom_interpolation.iteration(1,:);
            iteration_2=obj.fitted_geom_interpolation.iteration(2,:);
            
            index1 = permute(obj.fitted_mesh.connec(obj.cut_cells,iteration_1),[2 3 1]);
            index2 = permute(obj.fitted_mesh.connec(obj.cut_cells,iteration_2),[2 3 1]);
            gamma_1=obj.x_fitted(index1);
            gamma_2=obj.x_fitted(index2);
            coord1 = obj.fitted_mesh.coord(:,1); coord2 = obj.fitted_mesh.coord(:,2); coord3 = obj.fitted_mesh.coord(:,3);
            P1=[coord1(index1) coord2(index1) coord3(index1)];
            P2=[coord1(index2) coord2(index2) coord3(index2)];
            P=P1+gamma_1.*(P2-P1)./(gamma_1-gamma_2);
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

