classdef Filter_LevelSet_3D < Filter_LevelSet
    properties % For Facets plotting
        surrounding_facets_connectivities
    end
    
    methods
%         function setInterpolation_Unfitted(obj)
%             obj.interpolation_unfitted = Tetrahedra_Linear(obj.unfitted_mesh);
%         end
        
        function global_connectivities = computeFromLocalToGlobalConnectivities(obj,local_matrix_coordinates,global_matrix_coordinates,local_connectivities)
            indexes_in_global_matrix = obj.findIndexesOfCoordinatesAinCoordinateMatrixB(local_matrix_coordinates,global_matrix_coordinates);
            global_connectivities = indexes_in_global_matrix(local_connectivities);
        end
        
        function [boundary_facets_coordinates,boundary_facets_connectivities] = computeBoundaryFacets(obj,x)
            [interior_facets_coordinates, interior_facets_connectivities] = obj.computeInteriorFacets(x);
            [surrounding_active_facets_coordinates,surrounding_active_facets_connectivities] = obj.computeSurroundingActiveFacets(x);
            
            boundary_facets_coordinates = [surrounding_active_facets_coordinates;interior_facets_coordinates];
            boundary_facets_connectivities = [surrounding_active_facets_connectivities;interior_facets_connectivities+size(surrounding_active_facets_coordinates,1)];
        end
        
        function computeSurroundingFacets(obj)
            surrounding_facets_coordinates_raw = zeros(size(obj.mesh.coord)); surrounding_facets_connectivities_raw = zeros(size(obj.mesh.connec,1),obj.unfitted_mesh.ndim);
            k_coordinates = 0; k_connectivities = 0;
            for idime = 1:obj.unfitted_mesh.ndim
                [surrounding_facets_coordinates_raw, surrounding_facets_connectivities_raw,k_coordinates,k_connectivities] = obj.computeBoxFaceCoordNConnec(surrounding_facets_coordinates_raw,surrounding_facets_connectivities_raw,idime,max(obj.mesh.coord(:,idime)),k_coordinates,k_connectivities);
                [surrounding_facets_coordinates_raw, surrounding_facets_connectivities_raw,k_coordinates,k_connectivities] = obj.computeBoxFaceCoordNConnec(surrounding_facets_coordinates_raw,surrounding_facets_connectivities_raw,idime,min(obj.mesh.coord(:,idime)),k_coordinates,k_connectivities);
            end
            surrounding_facets_coordinates_raw(k_coordinates+1:end,:) = [];
            surrounding_facets_connectivities_raw(k_connectivities+1:end,:) = [];
            
            obj.surrounding_facets_connectivities = obj.computeFromLocalToGlobalConnectivities(surrounding_facets_coordinates_raw,obj.mesh.coord,surrounding_facets_connectivities_raw);
            
            %             figure
            %             patch('vertices',obj.mesh.coord,'faces',obj.surrounding_facets_connectivities,...
            %                 'edgecolor',[0.5 0 0], 'edgealpha',0.5,'edgelighting','flat',...
            %                 'facecolor','none','facelighting','flat')
            %             light
            %             axis equal off
        end
        
        function [surrounding_facets_coordinates, surrounding_facets_connectivities,k_coordinates,k_connectivities] = computeBoxFaceCoordNConnec(obj,surrounding_facets_coordinates,surrounding_facets_connectivities,idime,current_face_characteristic_coordinate,k_coordinates,k_connectivities)
            dimens = 1:obj.unfitted_mesh.ndim;
            new_coordinates = obj.mesh.coord(obj.mesh.coord(:,idime) == current_face_characteristic_coordinate,:);
            DT = delaunayTriangulation(new_coordinates(:,dimens(dimens ~= idime)));
            new_connectivities = DT.ConnectivityList + k_coordinates;
            surrounding_facets_coordinates(k_coordinates+1:k_coordinates+size(new_coordinates,1),:) = new_coordinates;
            surrounding_facets_connectivities(k_connectivities+1:k_connectivities+size(new_connectivities,1),:) = new_connectivities;
            k_coordinates = k_coordinates + size(new_coordinates,1);
            k_connectivities = k_connectivities + size(new_connectivities,1);
        end
        
        function [surrounding_active_facets_coordinates,surrounding_active_facets_connectivities] = computeSurroundingActiveFacets(obj,x)
            [full_elem,cut_elem]=obj.findCutElements(x,obj.surrounding_facets_connectivities);
            surrounding_active_full_facets_connectivities = obj.surrounding_facets_connectivities(full_elem,:);
            surrounding_active_full_facets_coordinates = obj.mesh.coord;
            
            
            
            %             hold on
            
            %             [interior_nodes_of_surrounding_active_cut_facets_coordinates,interior_nodes_of_surrounding_active_cut_facets_connectivities] = obj.computeSurfaceBoundaryInteriorSubcells(interp_element,obj.surrounding_facets_connectivities,x,cut_elem);
            
            %             plot3(interior_nodes_of_surrounding_active_cut_facets_coordinates(:,1),interior_nodes_of_surrounding_active_cut_facets_coordinates(:,2),interior_nodes_of_surrounding_active_cut_facets_coordinates(:,3),'.b')
            %
            %             surrounding_active_cut_facets_connectivities = obj.surrounding_facets_connectivities(cut_elem,:)
            
            interior_nodes_of_surrounding_active_cut_facets_coordinates = [];
            interior_nodes_of_surrounding_active_cut_facets_connectivities = [];
            
            surrounding_active_facets_coordinates = [surrounding_active_full_facets_coordinates;interior_nodes_of_surrounding_active_cut_facets_coordinates];
            surrounding_active_facets_connectivities = [surrounding_active_full_facets_connectivities;interior_nodes_of_surrounding_active_cut_facets_connectivities+size(surrounding_active_full_facets_coordinates,1)];
        end
        
        function [interior_facets_global_coordinates, interior_facets_global_connectivities] = computeInteriorFacets(obj,x)
            obj.unfitted_mesh = Mesh_Unfitted(obj.domainType,obj.mesh,obj.interpolation);        
            obj.unfitted_mesh.computeMesh(x);

            interior_facets_global_coordinates = obj.unfitted_mesh.coord;
            interior_facets_global_connectivities = obj.unfitted_mesh.connec;
        end
    end
    
    methods (Static, Access = public)                        
        function indexes = findIndexesOfCoordinatesAinCoordinateMatrixB(A,B)
            indexes = zeros(1,size(A,1));
            for inode = 1:size(A,1)
                match = true(size(B,1),1);
                for idime = 1:size(A,2)
                    match = match & B(:,idime) == A(inode,idime);
                end
                indexes(inode) = find(match,1);
            end
        end
    end
end

