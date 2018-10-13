classdef Filter_LevelSet_3D < Filter_LevelSet
    properties % For Facets plotting
        surrounding_facets_connectivities
    end
    
    methods
        function obj = Filter_LevelSet_3D
            obj.max_subcells = 20;
            obj.nnodes_subelem = 4;
            obj.ndim = 3;
        end
        
        function getQuadrature_Unfitted(obj)
            obj.quadrature_unfitted = Quadrature_Tetrahedra;
        end
        
        function getInterpolation_Unfitted(obj)
            obj.interpolation_unfitted = Tetrahedra_Linear(obj.unfitted_mesh);
        end
        
        function createUnfittedMesh_Interior(obj)
            obj.unfitted_mesh = Mesh_Unfitted_3D_Interior(obj.mesh.duplicate,obj.diffReacProb.geometry.interpolation);
        end
        
        function createUnfittedMesh_Boundary(obj)
            obj.unfitted_mesh = Mesh_Unfitted_3D_Boundary(obj.mesh.duplicate,obj.diffReacProb.geometry.interpolation);
        end
        
        function [interp_facet,quadrature_facet] = createFacet(obj)
            quadrature_facet = Quadrature.set('TRIANGLE');
            interp_facet = Triangle_Linear(obj.mesh);
            quadrature_facet.computeQuadrature(obj.quadrature_fitted.order);
            interp_facet.computeShapeDeriv(quadrature_facet.posgp);
        end
        
        function global_connectivities = computeFromLocalToGlobalConnectivities(obj,local_matrix_coordinates,global_matrix_coordinates,local_connectivities)
            indexes_in_global_matrix = obj.findCoordinatesIndexesInGlobalCoordinatesMatrix(local_matrix_coordinates,global_matrix_coordinates);
            global_connectivities = indexes_in_global_matrix(local_connectivities);
        end
        
        function [boundary_facets_coordinates,boundary_facets_connectivities] = computeBoundaryFacets(obj,x)
            [interior_facets_coordinates, interior_facets_connectivities] = obj.computeInteriorFacets(x);
            [surrounding_active_facets_coordinates,surrounding_active_facets_connectivities] = obj.computeSurroundingActiveFacets(x);
            
            boundary_facets_coordinates = [surrounding_active_facets_coordinates;interior_facets_coordinates];
            boundary_facets_connectivities = [surrounding_active_facets_connectivities;interior_facets_connectivities+size(surrounding_active_facets_coordinates,1)];
        end
        
        function computeSurroundingFacets(obj)
            surrounding_facets_coordinates_raw = zeros(size(obj.mesh.coord)); surrounding_facets_connectivities_raw = zeros(size(obj.mesh.connec,1),obj.ndim);
            k_coordinates = 0; k_connectivities = 0;
            for idime = 1:obj.ndim
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
            dimens = 1:obj.ndim;
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
            obj.createUnfittedMesh_Boundary;
            obj.unfitted_mesh.computeMesh(x);
            obj.unfitted_mesh.computeGlobalConnectivities;
            
            interior_facets_global_coordinates = obj.unfitted_mesh.coord;
            interior_facets_global_connectivities = obj.unfitted_mesh.connec;
        end
    end
    
    methods (Static)       
        function djacob = mapping(points,dvolu)
            v1 = diff(points([1 2],:));
            v2 = diff(points([1 3],:));
            A = 0.5*norm(cross(v1,v2));
            djacob = A/dvolu;
        end
        
        function indexes_in_global_matrix = findCoordinatesIndexesInGlobalCoordinatesMatrix(coordinates_local,coordinates_global)
            indexes_in_global_matrix = zeros(1,size(coordinates_local,1));
            for inode = 1:size(coordinates_local,1)
                match = true(size(coordinates_global,1),1);
                for idime = 1:size(coordinates_local,2)
                    match = match & coordinates_global(:,idime) == coordinates_local(inode,idime);
                end
                indexes_in_global_matrix(inode) = find(match,1);
            end
        end
    end
end

