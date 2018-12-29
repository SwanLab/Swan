classdef CutPoints_Calculator_3D < CutPoints_Calculator_Abstract
    methods (Access = public, Static)
        function [P,active_nodes] = computeCutPoints_Iso(mesh_background,x_background,background_cut_cells,background_geom_interpolation)
            pos_nodes = background_geom_interpolation.pos_nodes;
            
            iteration_1 = background_geom_interpolation.iteration(1,:);
            iteration_2 = background_geom_interpolation.iteration(2,:);
            
            gamma_1 = permute(x_background(mesh_background.connec(background_cut_cells,iteration_1)),[2 3 1]);
            gamma_2 = permute(x_background(mesh_background.connec(background_cut_cells,iteration_2)),[2 3 1]);
            
            P1 = repmat(pos_nodes(iteration_1,:),[1 1 size(background_cut_cells)]);
            P2 = repmat(pos_nodes(iteration_2,:),[1 1 size(background_cut_cells)]);
            P = P1+gamma_1.*(P2-P1)./(gamma_1-gamma_2);
            
            active_nodes = sign(gamma_1.*gamma_2)<=0;
        end
        
        function [P,active_nodes] = computeCutPoints_Global(mesh_background,x_background,background_cut_cells,background_geom_interpolation)
            iteration_1 = background_geom_interpolation.iteration(1,:);
            iteration_2 = background_geom_interpolation.iteration(2,:);
            
            index1 = permute(mesh_background.connec(background_cut_cells,iteration_1),[2 3 1]);
            index2 = permute(mesh_background.connec(background_cut_cells,iteration_2),[2 3 1]);
            
            gamma_1 = x_background(index1);
            gamma_2 = x_background(index2);
            
            coord1 = mesh_background.coord(:,1);
            coord2 = mesh_background.coord(:,2);
            coord3 = mesh_background.coord(:,3);
            
            P1 = [coord1(index1) coord2(index1) coord3(index1)];
            P2 = [coord1(index2) coord2(index2) coord3(index2)];
            P = P1+gamma_1.*(P2-P1)./(gamma_1-gamma_2);
            
            active_nodes = sign(gamma_1.*gamma_2)<=0;
        end
    end
end

