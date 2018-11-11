classdef Mesh_Unfitted_2D_Delaunay < Mesh_Unfitted_2D & Mesh_Unfitted_Delaunay
    
    methods
        function obj = Mesh_Unfitted_2D_Delaunay(fitted_mesh,x_fitted,fitted_geom_interpolation)
            obj@Mesh_Unfitted_2D(fitted_mesh,x_fitted,fitted_geom_interpolation);
            obj.max_subcells = 6;
            obj.nnodes_subcell = 3;
        end        
        function [P,active_nodes] = findCutPoints_Iso(obj)
            pos_nodes = obj.fitted_geom_interpolation.pos_nodes;
            
            gamma_1 = permute(obj.x_fitted(obj.fitted_mesh.connec(obj.cut_cells,:)),[2 3 1]);
            gamma_2 = permute([obj.x_fitted(obj.fitted_mesh.connec(obj.cut_cells,2:end)),obj.x_fitted(obj.fitted_mesh.connec(obj.cut_cells,1))],[2 3 1]);
            
            P1 = repmat(pos_nodes,[1 1 size(obj.cut_cells)]);
            % !! There's an specific Matlab function to do this patch!!
            P2 = repmat([pos_nodes(2:end,:);pos_nodes(1,:)],[1 1 size(obj.cut_cells)]);
            P = P1+gamma_1.*(P2-P1)./(gamma_1-gamma_2);
            
            active_nodes = sign(gamma_1.*gamma_2)<=0;
        end       
    end
end

