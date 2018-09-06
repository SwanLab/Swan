classdef Mesh_Unfitted_3D_Facets < Mesh_Unfitted_3D
    methods
        function obj = Mesh_Unfitted_3D_Facets(fitted_mesh,x_fitted,fitted_geom_interpolation)
            obj@Mesh_Unfitted_3D(fitted_mesh,x_fitted,fitted_geom_interpolation);
        end
        
        function boundary_subfacets_connectivities = computeFacetsLocalConnectivities(obj,cutPoints_iso,pos_nodes, phi)
            del_coord = [pos_nodes; cutPoints_iso];
            DT = delaunayTriangulation(del_coord);
            subcells_connectivities = DT.ConnectivityList;
            
            interior_nodes = find(phi<=0); exterior_nodes = find(phi>0);
            
            % Find subcells formed by 1 interior node & 3 cutPoints
            counter_interior_nodes = obj.countInputNodesPerCell(subcells_connectivities,interior_nodes); %#ok<FNDSB>
            counter_exterior_nodes = obj.countInputNodesPerCell(subcells_connectivities,exterior_nodes); %#ok<FNDSB>
            boundary_subcells_connectivities = subcells_connectivities(counter_interior_nodes == 1 & counter_exterior_nodes == 0,:);
            
            boundary_subfacets_connectivities = zeros([size(boundary_subcells_connectivities,1),3]);
            for i = 1:size(boundary_subcells_connectivities,1)
                boundary_subfacets_connectivities(i,:) = boundary_subcells_connectivities(i,boundary_subcells_connectivities(i,:)>size(pos_nodes,1));
            end
            boundary_subfacets_connectivities = boundary_subfacets_connectivities - size(pos_nodes,1);
        end
        
        function global_connectivities = computeFromLocalToGlobalConnectivities(obj,local_matrix_coordinates,global_matrix_coordinates,local_connectivities)
            indexes_in_global_matrix = obj.findCoordinatesIndexesInGlobalCoordinatesMatrix(local_matrix_coordinates,global_matrix_coordinates);
            global_connectivities = indexes_in_global_matrix(local_connectivities);
        end
    end
    
    methods (Static)
        function all_cutPoints_global = findActiveCutPoints(P_global,active_nodes_global)
            P_global_x = P_global(:,1,:); P_global_y = P_global(:,2,:); P_global_z = P_global(:,3,:);
            all_cutPoints_global(:,1) = P_global_x(active_nodes_global); all_cutPoints_global(:,2) = P_global_y(active_nodes_global); all_cutPoints_global(:,3) = P_global_z(active_nodes_global);
            
            all_cutPoints_global = unique(all_cutPoints_global,'rows','stable');
        end
        
        function counter = countInputNodesPerCell(connectivities,nodes)
            counter = zeros(size(connectivities,1),1);
            for iconnec = 1:size(nodes,1)
                match = false(size(connectivities,1),1);
                for inode = 1:size(connectivities,2)
                    match = match | connectivities(:,inode) == nodes(iconnec);
                end
                counter = counter + match ;
            end
        end
    end
end

