classdef Mesh_Unfitted_3D_Boundary < Mesh_Unfitted_3D & Mesh_Unfitted_Boundary
    methods
        function obj = Mesh_Unfitted_3D_Boundary(fitted_mesh,fitted_geom_interpolation)
            obj.storeFittedMesh(fitted_mesh,fitted_geom_interpolation);
            obj.max_subcells = 6; % !! ?? !!
            obj.nnodes_subcell = 3;
        end
        
        %         function obj = computeMesh_withExtBounds(obj,x)
        %             interior_boundary_mesh = Mesh_Unfitted_3D_Boundary(obj.fitted_mesh.duplicate,obj.fitted_geom_interpolation);
        %             interior_boundary_mesh.computeMesh(x);
        %             surrounding_boundary_meshes = obj.computeSurrondingBoundaryMeshes;
        %         end
        
        function plot(obj)
            figure, hold on
            patch('vertices',obj.coord,'faces',obj.connec,...
                'edgecolor',[0.5 0 0], 'edgealpha',0.5,'edgelighting','flat',...
                'facecolor',[1 0 0],'facelighting','flat')
            light
            axis equal off
            hold off
        end
    end
    
    methods (Access = ?Mesh_Unfitted_Boundary)
        function facets_connec = computeFacetsConnectivities(obj,~,interior_subcell_coord_iso,cell_x_value)
            subcells_connec = obj.computeDelaunay(interior_subcell_coord_iso);
            boundary_subcells_connec = obj.findBoundarySubcells(subcells_connec,cell_x_value);
            
            number_nodes = size(obj.fitted_geom_interpolation.pos_nodes,1);
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
        %                     surrounding_boundary_meshes{iface} =  Mesh_Unfitted_2D_Interior(obj.fitted_mesh.duplicate,geom_interpolation);
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
end
