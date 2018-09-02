classdef Mesh_Unfitted_2D < Mesh_Unfitted
    
    methods
        function obj = Mesh_Unfitted_2D(fitted_mesh,x_fitted,fitted_geom_interpolation)
            obj@Mesh_Unfitted(fitted_mesh,fitted_geom_interpolation,x_fitted);
            obj.max_subcells = 6;
            obj.nnodes_subcell = 3;
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
        
        function computeDvoluCut(obj)
            x1 = obj.coord_cut(:,1,1); y1 = obj.coord_cut(:,1,2);
            x2 = obj.coord_cut(:,2,1); y2 = obj.coord_cut(:,2,2);
            x3 = obj.coord_cut(:,3,1); y3 = obj.coord_cut(:,3,2);
            
            obj.dvolu_cut= 0.5*abs((x2-x1).*(y3-y1)-(x3-x1).*(y2-y1));
        end
    end
end

