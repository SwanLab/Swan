classdef Mesh_Unfitted_Composite < Mesh_Unfitted_Abstract
    properties (Access = private)
        meshBackground;
        meshInterior
        boxFaceMeshes
        nodesInBoxFaces
        activeBoxFaceMesh
        
        ndim
        nsides = 2;
        nboxFaces
    end
    
    methods (Access = public)
        function obj = Mesh_Unfitted_Composite(meshType,meshBackground,interpolation_background)
            obj.init(meshBackground);
            obj.createInteriorMesh(meshType,meshBackground,interpolation_background)
            obj.createBoxMeshes()
        end
        
        function computeMesh(obj,levelSet)
            obj.computeInteriorMesh(levelSet);
            obj.computeBoxMeshes(levelSet);
        end
        
        function M = computeMass(obj)
            Mi = obj.computeInteriorMass();
            Mb = obj.computeBoxMass();
            M = Mi + Mb;
        end
    end
    
    methods (Access = private)
        function init(obj,meshBackground)
            obj.ndim = meshBackground.ndim;
            obj.meshBackground = meshBackground;
            obj.nboxFaces = obj.ndim*obj.nsides;
        end
        
        function createInteriorMesh(obj,meshType,meshBackground,interpolation_background)
            obj.meshInterior = Mesh_Unfitted(meshType,meshBackground,interpolation_background);
        end
        
        function computeInteriorMesh(obj,levelSet)
            obj.meshInterior.computeMesh(levelSet);
        end
        
        function M = computeInteriorMass(obj)
            M = obj.meshInterior.computeMass();
        end
        
        function createBoxMeshes(obj)
            iface = 0;
            for idime = 1:obj.ndim
                for iside = 1:obj.nsides
                    iface = iface + 1;
                    [boxFaceMesh,nodesInBoxFace] = obj.createBoxFaceMesh(idime,iside);
                    obj.boxFaceMeshes{iface}     = boxFaceMesh;
                    obj.nodesInBoxFaces{iface}   = nodesInBoxFace;
                end
            end
        end
        
        function computeBoxMeshes(obj,levelSet)
            iface = 0;
            obj.activeBoxFaceMesh = false([1 obj.nboxFaces]);
            for idime = 1:obj.ndim
                for iside = 1:obj.nsides
                    iface = iface + 1;
                    boxFaceMesh = obj.boxFaceMeshes{iface};
                    mshBack = boxFaceMesh.meshBackground;
                    lsBoxFace = levelSet(obj.nodesInBoxFaces{iface});
                    if obj.isBoxMeshActive(lsBoxFace,mshBack)
                        obj.boxFaceMeshes{iface}.computeMesh(lsBoxFace);
                        obj.activeBoxFaceMesh(iface) = true;
                    end
                end
            end
        end
        
        function M = computeBoxMass(obj)
            M = 0;
            for iface = 1:obj.nboxFaces
                if obj.activeBoxFaceMesh(iface)
                    M = M + obj.boxFaceMeshes{iface}.computeMass();
                end
            end
        end
        
        function [boxFaceMesh,nodesInBoxFace] = createBoxFaceMesh(obj,idime,iside)
            [mb,nodesInBoxFace] = obj.createBoxFaceBackgroundMesh(idime,iside);
            interpolation_unfitted = Interpolation.create(mb,'LINEAR');
            boxFaceMesh = Mesh_Unfitted('INTERIOR',mb,interpolation_unfitted);
        end
        
        function [mb, nodesInBoxFace] = createBoxFaceBackgroundMesh(obj,idime,iside)
            [boxFaceCoords,nodesInBoxFace] = obj.getFaceCoordinates(idime,iside);
            DT = delaunayTriangulation(boxFaceCoords);
            face_connec = DT.ConnectivityList;
            mb = Mesh;
            mb = mb.create(boxFaceCoords,face_connec);
        end
        
        function [boxFaceCoords, nodesInBoxFace] = getFaceCoordinates(obj,idime,iside)
            if iside == 1
                L = min(obj.meshBackground.coord(:,idime));
            else
                L = max(obj.meshBackground.coord(:,idime));
            end
            
            nodesInBoxFace = obj.meshBackground.coord(:,idime) == L;
            boxFaceCoords = obj.meshBackground.coord(nodesInBoxFace,:);
            boxFaceCoords = obj.removeExtraDimension(boxFaceCoords,idime);
        end
        
        function indexes = findConnecIndexes(obj,coord_indexes,nnode)
            number_of_valid_nodes_per_element = sum(ismember(obj.meshBackground.connec,coord_indexes),2);
            indexes = number_of_valid_nodes_per_element == nnode;
        end
    end
    
    methods (Static, Access = private)
        function itIs = isBoxMeshActive(levelSet,meshBack)
            phi_nodes = levelSet(meshBack.connec);
            phi_case = sum((sign(phi_nodes)<0),2);
            itIs = (any(phi_case));
        end
        
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

