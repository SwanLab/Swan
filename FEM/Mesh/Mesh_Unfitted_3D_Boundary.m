classdef Mesh_Unfitted_3D_Boundary < Mesh_Unfitted_3D & Mesh_Unfitted_Boundary
    methods
        function obj = Mesh_Unfitted_3D_Boundary(fitted_mesh,fitted_geom_interpolation)
            obj.storeFittedMesh(fitted_mesh,fitted_geom_interpolation);
            obj.geometryType = 'TRIANGLE';
            obj.max_subcells = 6; % !! ?? !!
            obj.nnodes_subcell = 3;
        end
        
        function [facets_coord_iso,facets_coord_global,facets_x_value,facets_connec] = computeSubcells(obj,fitted_cell_connec,subcell_cutPoints_iso,subcell_cutPoints_global)
            subcell_coord_iso = [obj.fitted_geom_interpolation.pos_nodes; subcell_cutPoints_iso];
%             subcell_coord_global = [obj.fitted_mesh.coord(fitted_cell_connec,:); subcell_cutPoints_global];
            subcell_x_value = obj.x_fitted(fitted_cell_connec)'; % !! RENAME !!
            
            % !! MOVE THIS TO THE BOUNDARY SUPERCLASS !!
            facets_coord_iso = subcell_cutPoints_iso;
            facets_coord_global = subcell_cutPoints_global;
            facets_x_value = zeros(1,size(subcell_cutPoints_iso,1)); % !! RENAME !!
            
            facets_connec = obj.computeFacetsConnectivities(subcell_coord_iso,subcell_x_value);
        end
        
        function plot(obj)
            figure, hold on
            patch('vertices',obj.unfitted_coord_global,'faces',obj.unfitted_connec_global,...
                'edgecolor',[0.5 0 0], 'edgealpha',0.5,'edgelighting','flat',...
                'facecolor',[1 0 0],'facelighting','flat')
            light
            axis equal off
            hold off
        end
    end
    
    methods (Access = private)
        function facets_connec = computeFacetsConnectivities(obj,subcell_coord_iso,phi)
            number_nodes = size(obj.fitted_geom_interpolation.pos_nodes,1);
            DT = delaunayTriangulation(subcell_coord_iso);
            subcells_connec = DT.ConnectivityList;
            
            interior_nodes = find(phi<=0); exterior_nodes = find(phi>0);
            
            % Find subcells formed by 1 interior node & 3 cutPoints (PACK IN FUNCTION NAMED LIKE THIS)
            number_interior_nodes = obj.countInputNodesPerCell(subcells_connec,interior_nodes); %#ok<FNDSB>
            number_exterior_nodes = obj.countInputNodesPerCell(subcells_connec,exterior_nodes); %#ok<FNDSB>
            boundary_subcells_connec = subcells_connec(number_interior_nodes == 1 & number_exterior_nodes == 0,:);
            
            facets_connec = zeros([size(boundary_subcells_connec,1),3]);
            for i = 1:size(boundary_subcells_connec,1)
                facets_connec(i,:) = boundary_subcells_connec(i,boundary_subcells_connec(i,:)>number_nodes);
            end
            
            % !!!!!!!!!!!!!!!!!!!!!!! PLOTTING !!!!!!!!!!!!!!!!!!!!!!!!
            %             figure, hold on
            %             fac = [1 2 3 4; 2 6 7 3; 4 3 7 8; 1 5 8 4; 1 2 6 5; 5 6 7 8];
            %             patch('Faces',fac,'Vertices',[-1 -1 -1; -1 1 -1; 1 1 -1; 1 -1 -1; -1 -1 1; -1 1 1; 1 1 1; 1 -1 1],'FaceColor','w','FaceAlpha',0.0);
            %             for iconnec = 1:size(boundary_subfacets_connec,1)
            %                 plot3(delaunay_coord(boundary_subfacets_connec(iconnec,:),1),delaunay_coord(boundary_subfacets_connec(iconnec,:),2),delaunay_coord(boundary_subfacets_connec(iconnec,:),3),'.-g')
            %                 plot3(delaunay_coord(boundary_subfacets_connec(iconnec,[1 3]),1),delaunay_coord(boundary_subfacets_connec(iconnec,[1 3]),2),delaunay_coord(boundary_subfacets_connec(iconnec,[1 3]),3),'.-g')
            %             end
            %             patch('vertices',delaunay_coord,'faces',boundary_subfacets_connec,'edgecolor',[0 0.5 0],...
            %                 'facecolor',[0 1 0],'facelighting','phong')
            %
            %             plot3(cutPoints_iso(:,1),cutPoints_iso(:,2),cutPoints_iso(:,3),'og')
            %
            %             plot3(pos_nodes(phi>0,1),pos_nodes(phi>0,2),pos_nodes(phi>0,3),'+b')
            %             plot3(pos_nodes(phi<0,1),pos_nodes(phi<0,2),pos_nodes(phi<0,3),'<b')
            %
            %             axis equal
            %             view([115 20])
            %
            %             close
            facets_connec = facets_connec - number_nodes;
        end
    end
    
    methods (Access = private, Static)
        function counter = countInputNodesPerCell(connectivities,nodes)
            counter = zeros(size(connectivities,1),1);
            for inode = 1:length(nodes)
                match = false(size(connectivities,1),1);
                for iconnec = 1:size(connectivities,2)
                    match = match | connectivities(:,iconnec) == nodes(inode);
                end
                counter = counter + match ;
            end
        end
    end
end
