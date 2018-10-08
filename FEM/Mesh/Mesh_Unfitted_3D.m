classdef Mesh_Unfitted_3D < Mesh_Unfitted
    methods
        function obj = Mesh_Unfitted_3D(fitted_mesh,x_fitted,fitted_geom_interpolation)
            obj@Mesh_Unfitted(fitted_mesh,x_fitted,fitted_geom_interpolation);
            obj.max_subcells = 20;
            obj.nnodes_subcell = 4;
        end
        
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
        
        function computeDvoluCut(obj)
            x1 = obj.unfitted_cut_coord_iso_per_cell(:,1,1);   y1 = obj.unfitted_cut_coord_iso_per_cell(:,1,2);   z1= obj.unfitted_cut_coord_iso_per_cell(:,1,3);
            x2 = obj.unfitted_cut_coord_iso_per_cell(:,2,1);   y2 = obj.unfitted_cut_coord_iso_per_cell(:,2,2);   z2= obj.unfitted_cut_coord_iso_per_cell(:,2,3);
            x3 = obj.unfitted_cut_coord_iso_per_cell(:,3,1);   y3 = obj.unfitted_cut_coord_iso_per_cell(:,3,2);   z3= obj.unfitted_cut_coord_iso_per_cell(:,3,3);
            x4 = obj.unfitted_cut_coord_iso_per_cell(:,4,1);   y4 = obj.unfitted_cut_coord_iso_per_cell(:,4,2);   z4= obj.unfitted_cut_coord_iso_per_cell(:,4,3);
            
            J = x1.*y3.*z2-x1.*y2.*z3+x2.*y1.*z3-x2.*y3.*z1-x3.*y1.*z2+x3.*y2.*z1+x1.*y2.*z4-x1.*y4.*z2-x2.*y1.*z4+x2.*y4.*z1+...
                x4.*y1.*z2-x4.*y2.*z1-x1.*y3.*z4+x1.*y4.*z3+x3.*y1.*z4-x3.*y4.*z1-x4.*y1.*z3+x4.*y3.*z1+x2.*y3.*z4-x2.*y4.*z3...
                -x3.*y2.*z4+x3.*y4.*z2+x4.*y2.*z3-x4.*y3.*z2;
            obj.dvolu_cut = J/6;
        end
        
        function assignUnfittedCutCoordIsoPerCell(obj,new_unfitted_cut_coord_iso,new_subcell_cut_interior_connec_iso,c0,c1)
            new_unfitted_cut_coord_iso_x = new_unfitted_cut_coord_iso(:,1);
            new_unfitted_cut_coord_iso_y = new_unfitted_cut_coord_iso(:,2);
            new_unfitted_cut_coord_iso_z = new_unfitted_cut_coord_iso(:,3);
            
            obj.unfitted_cut_coord_iso_per_cell(c0+1:c1,:,1) = new_unfitted_cut_coord_iso_x(new_subcell_cut_interior_connec_iso);
            obj.unfitted_cut_coord_iso_per_cell(c0+1:c1,:,2) = new_unfitted_cut_coord_iso_y(new_subcell_cut_interior_connec_iso);
            obj.unfitted_cut_coord_iso_per_cell(c0+1:c1,:,3) = new_unfitted_cut_coord_iso_z(new_subcell_cut_interior_connec_iso);
        end
    end
end

