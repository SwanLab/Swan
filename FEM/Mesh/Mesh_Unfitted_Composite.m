classdef Mesh_Unfitted_Composite < Mesh
    properties (GetAccess = public, SetAccess = private)
        interior_mesh
        surrounding_meshes
    end
    
    methods (Access = public, Static)
        function obj = create(mesh_background,interpolation_background)
            obj = Mesh_Unfitted_Composite;
            obj.interior_mesh = Mesh_Unfitted_Factory.create('INTERIOR',mesh_background,interpolation_background);
        end
    end
    
    methods (Access = public)
        function computeMesh(obj,x)
            obj.interior_mesh.computeMesh(x);
            obj.computeSurroundingMeshes(x);
        end
    end
    
    methods (Access = private)
        function computeSurroundingMeshes(obj)
            interpolation_unfitted = Interpolation.create(obj.interior_mesh,'LINEAR');
            for idime = 1:obj.ndim
                for iside = 1:2
                    [face_mesh,valid_nodes] = obj.createFaceMesh(idime,iside);
                    unfitted_mesh2D = Mesh_Unfitted_2D_Interior(face_mesh,interpolation_unfitted);
                    unfitted_mesh2D.computeMesh(obj.x_background(valid_nodes));
                    
                    if ~isempty(unfitted_mesh2D.connec)
                        
                    end
                end
            end
        end
        
        function [face_mesh, valid_nodes] = createFaceMesh(obj,idime,iside)
            [face_coord,valid_nodes] = obj.getFaceCoordinates(idime,iside);
            % face_connec = obj.getFaceConnectivities(face_coord,idime,iside);
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
            % face_coord = unique(face_coord,'row');
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

