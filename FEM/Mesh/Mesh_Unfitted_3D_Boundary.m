classdef Mesh_Unfitted_3D_Boundary < Mesh_Unfitted_3D & Mesh_Unfitted_Boundary
    methods (Access = public)
        function obj = Mesh_Unfitted_3D_Boundary(mesh_background,background_geom_interpolation)
            obj.storeBackgroundMesh(mesh_background,background_geom_interpolation);
            obj.max_subcells = 6; % !! ?? !!
            obj.nnodes_subcell = 3;
        end
        
        %         function obj = computeMesh_withExtBounds(obj,x)
        %             interior_boundary_mesh = Mesh_Unfitted_3D_Boundary(obj.mesh_background.clone,obj.background_geom_interpolation);
        %             interior_boundary_mesh.computeMesh(x);
        %             surrounding_boundary_meshes = obj.computeSurrondingBoundaryMeshes;
        %         end
        
        function add2plot(obj,h)
            hold on;
            patch(h,'vertices',obj.coord,'faces',obj.connec,...
                'edgecolor',[0.5 0 0], 'edgealpha',0.5,'edgelighting','flat',...
                'facecolor',[1 0 0],'facelighting','flat')
        end
        
        function S = computeSurface(obj)
            S = obj.computeMass;
            
            interpolation_unfitted = Interpolation.create(obj,'LINEAR');
            for idime = 1:obj.ndim
                for iside = 1:2
                    [face_mesh,valid_nodes] = obj.createFaceMesh(idime,iside);
                    unfitted_mesh2D = Mesh_Unfitted_2D_Interior(face_mesh,interpolation_unfitted);
                    unfitted_mesh2D.computeMesh(obj.x_background(valid_nodes));
                    
                    if ~isempty(unfitted_mesh2D.connec)
                        S_2D = unfitted_mesh2D.computeSurface;
                        S = S + S_2D;
                    end
                end
            end
        end
    end
    
    methods (Access = ?Mesh_Unfitted_Boundary)
        function facets_connec = computeFacetsConnectivities(obj,~,interior_subcell_coord_iso,cell_x_value)
            subcells_connec = obj.computeDelaunay(interior_subcell_coord_iso);
            boundary_subcells_connec = obj.findBoundarySubcells(subcells_connec,cell_x_value);
            
            number_nodes = size(obj.background_geom_interpolation.pos_nodes,1);
            facets_connec = zeros([size(boundary_subcells_connec,1),3]);
            for i = 1:size(boundary_subcells_connec,1)
                facets_connec(i,:) = boundary_subcells_connec(i,boundary_subcells_connec(i,:)>number_nodes);
            end
            facets_connec = facets_connec - number_nodes;
        end
    end
    
    methods (Access = private)
        function boundary_subcells_connec = findBoundarySubcells(obj,interior_subcells_connec,phi)
            % Find subcells formed by 1 interior node & 3 cutPoints (0 exterior nodes)
            interior_nodes = find(phi<=0);
            exterior_nodes = find(phi>0);
            
            number_interior_nodes = obj.countInputNodesPerCell(interior_subcells_connec,interior_nodes); %#ok<FNDSB>
            number_exterior_nodes = obj.countInputNodesPerCell(interior_subcells_connec,exterior_nodes); %#ok<FNDSB>
            boundary_subcells_connec = interior_subcells_connec(number_interior_nodes == 1 & number_exterior_nodes == 0,:);
        end
        
        %         function surrounding_boundary_meshes = computeSurrondingBoundaryMeshes(obj)
        % %             geom_interpolation = Geometry(..);
        %             nfaces = 2*obj.ndim;
        %             surrounding_boundary_meshes = cell([1 nfaces]);
        %             domain_limits = obj.getDomainLimits;
        %             for idime = 1:obj.ndim
        %                 for iside = 1:2
        %                     iface = 2*(idime-1) + iside;
        %                     face_coord = obj.getFaceCoordinates(domain_limits(idime,iside),idime);
        %                     surrounding_boundary_meshes{iface} =  Mesh_Unfitted_2D_Interior(obj.mesh_background.clone,geom_interpolation);
        %                 end
        %             end
        %         end
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
    
    %% !! DIRECTLY MOVED FROM FILTER_LEVELSET_3D_Boundary !!
    methods ( Access = private)
        function [face_mesh, valid_nodes] = createFaceMesh(obj,idime,iside)
            [face_coord,valid_nodes] = obj.getFaceCoordinates(idime,iside);
            %             face_connec = obj.getFaceConnectivities(face_coord,idime,iside);
            face_connec = obj.computeDelaunay(face_coord);
            face_mesh = Mesh;
            face_mesh = face_mesh.create(face_coord,face_connec);
        end
        
        function [face_coord, valid_coord] = getFaceCoordinates(obj,idime,iside)
            if iside == 1
                L = min(obj.mesh_background.coord(:,idime));
            else
                L = max(obj.mesh_background.coord(:,idime));
            end
            
            valid_coord = obj.mesh_background.coord(:,idime) == L;
            face_coord = obj.mesh_background.coord(valid_coord,:);
            face_coord = obj.removeExtraDimension(face_coord,idime);
            %             face_coord = unique(face_coord,'row');
        end
        
        %         function face_connec = getFaceConnectivities(obj,face_coord,idime,iside)
        %             indexes_in_global_matrix = obj.findIndexesOfCoordinatesAinCoordinateMatrixB(face_coord,obj.mesh_background.coord);
        %             nnode = 3;
        %             valid_cells = sum(ismember(obj.mesh_background.connec,indexes_in_global_matrix),2)==nnode;
        %             face_connec = obj.mesh_background.connec(valid_cells,:);
        %             %!! NOT WORKING !!
        %         end
        
        function indexes = findConnecIndexes(obj,coord_indexes,nnode)
            number_of_valid_nodes_per_element = sum(ismember(obj.mesh_background.connec,coord_indexes),2);
            indexes = number_of_valid_nodes_per_element == nnode;
        end
    end
    
    methods (Static, Access = private)
        function face_connec = removeExtraNodes(face_connec_raw,coord_indexes,nnode)
            valid_nodes = ismember(face_connec_raw,coord_indexes);
            
            face_connec = zeros(size(face_connec_raw,1),nnode);
            for i = 1:size(face_connec,1)
                face_connec(i,:) = face_connec_raw(i,valid_nodes(i,:));
            end
        end
        
        function face_coord = removeExtraDimension(face_coord,idime)
            dimen = [1 2 3];
            face_coord = face_coord(:,dimen(dimen~=idime));
        end
    end
end
