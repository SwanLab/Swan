classdef Mesh_Unfitted_2D_Boundary < Mesh_Unfitted_2D & Mesh_Unfitted_Boundary
    methods
        function obj = Mesh_Unfitted_2D_Boundary(fitted_mesh,fitted_geom_interpolation)
            obj.storeFittedMesh(fitted_mesh,fitted_geom_interpolation);
            obj.geometryType = 'LINE';
            obj.max_subcells = 2;
            obj.nnodes_subcell = 2;
        end
        
        function [subcell_coord_iso,subcell_coord_global,subcell_x_value,interior_subcell_connec] = computeSubcells(obj,~,subcell_cutPoints_iso,subcell_cutPoints_global)
            subcell_coord_iso = subcell_cutPoints_iso;
            subcell_coord_global = subcell_cutPoints_global;
            subcell_x_value = zeros(1,size(subcell_cutPoints_iso,1));
            
            interior_subcell_connec = obj.computeBoundarySubcellsConnectivities(subcell_coord_iso);
        end
        
        function plot(obj)
            figure, hold on
            for icell = 1:size(obj.unfitted_connec_global,1)
                plot(obj.unfitted_coord_global(obj.unfitted_connec_global(icell,:),1),obj.unfitted_coord_global(obj.unfitted_connec_global(icell,:),2),'k-');
            end
            axis equal off
            hold off
        end
    end
    
    methods (Static, Access = private)
        function subcell_connec = computeBoundarySubcellsConnectivities(subcell_coord_iso)
            if size(subcell_coord_iso,1) == 2
                subcell_connec = [1 2];
            elseif size(subcell_coord_iso,1) == 4
                DT = delaunayTriangulation(subcell_coord_iso);
                del_connec = DT.ConnectivityList;
                
                node_positive_iso = find(obj.x_fitted(inode_global)>0);
                %                 subcell_connec = zeros(length(node_positive_iso),size(del_connec,2));
                for idel = 1:length(node_positive_iso)
                    [connec_positive_nodes, ~] = find(del_connec==node_positive_iso(idel));
                    subcell_connec(idel,:) = del_connec(connec_positive_nodes(end),del_connec(connec_positive_nodes(end),:)~=node_positive_iso(idel))-interpolation.nnode;
                end
            else
                error('Case not considered.')
            end
        end
    end
end