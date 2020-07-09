classdef Mesh_Total < Mesh_Composite 
    
    properties (GetAccess = public, SetAccess = private)
        nBoxFaces
        
        innerMeshOLD
        boxFaceMeshes
        nodesInBoxFaces
        
        removedDimensions
        removedDimensionCoord
        
        npnod
        nnode
        embeddedDim

    end
    
    properties (Access = private)
        borderNodes
        borderElements
        isExteriorMeshExplicit
    end
    
    methods (Access = public)
        
        function obj = Mesh_Total(cParams)
            obj.init(cParams);
            obj.createInteriorMesh();
            obj.createBoxFaceMeshes();
            obj.defineActiveMeshes();            
            obj.geometryType = obj.innerMeshOLD.geometryType;
            obj.nelem = size(obj.connec,1);
            obj.npnod = obj.innerMeshOLD.npnod;
            obj.nnode = obj.innerMeshOLD.nnode;
            obj.createInterpolation();
            obj.computeElementCoordinates();             
        end
        
        function S = computeMeanCellSize(obj)
            S = obj.innerMeshOLD.computeMeanCellSize();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.coord  = cParams.coord;
            obj.connec = cParams.connec;
            obj.obtainExteriorMesh(cParams);              
            obj.ndim   = size(obj.coord,2);
            obj.embeddedDim = obj.ndim;
        end
        
        function obtainExteriorMesh(obj,cParams)
            obj.isExteriorMeshExplicit = false;
            if isfield(cParams,'borderNodes')
                if ~isempty(cParams.borderNodes)
                    obj.borderNodes    = cParams.borderNodes;
                    obj.borderElements = cParams.borderElements;
                    obj.isExteriorMeshExplicit = true;
                end
            end
        end
        
        function createInteriorMesh(obj)
            s.connec = obj.connec;
            s.coord  = obj.coord;
            obj.innerMeshOLD = Mesh().create(s);
            obj.append(obj.innerMeshOLD);
        end
        
        function createBoxFaceMeshes(obj)
            if obj.isExteriorMeshExplicit
                obj.computeExteriorMeshesFromData();
            else
               obj.computeExteriorMeshesFromBoxSides();
            end
        end
        
        function computeExteriorMeshesFromData(obj)
            nExteriorMeshes = 1;
            for imesh = 1:nExteriorMeshes
                nodes  = obj.borderNodes;
                s.coord = obj.coord(nodes,:);
                s.connec = obj.computeConnectivitiesFromData();
                s.type = 'BOUNDARY';
                m = Mesh().create(s);
                obj.boxFaceMeshes{imesh} = m;
                obj.append(m);
                obj.nodesInBoxFaces{imesh} = false(size(obj.coord,1),1);
                obj.nodesInBoxFaces{imesh}(nodes,1) = true;
                obj.globalConnectivities{imesh} = obj.borderElements(:,2:end);  
            end
           obj.nBoxFaces = numel(obj.boxFaceMeshes);                              
        end
        
        function borderConnecSwitch = computeConnectivitiesFromData(obj)
            connec = obj.borderElements(:,2:end);            
            nNode = size(connec,2);
            nElem = size(connec,1);
            icell = 1;
            borderConnec(1,:) = connec(1,:);
            nodeOld = connec(icell,1);
            nodeNew = connec(icell,2);            
            for ielem = 2:nElem 
                isInElementAsNodeA = find(connec(:,1) == nodeNew);
                isInElementAsNodeB = find(connec(:,2) == nodeNew);
                
                isNewCellAsA = isInElementAsNodeA ~= icell;
                isNewCellAsB = isInElementAsNodeB ~= icell;
                
                if isNewCellAsA                    
                    icell = isInElementAsNodeA;                    
                    nodeOld = nodeNew;                    
                    nodeNew = connec(icell,2);
                    borderConnec(ielem,1) = nodeOld;
                    borderConnec(ielem,2) = nodeNew;
                elseif isNewCellAsB
                    icell = isInElementAsNodeB;                    
                    nodeOld = nodeNew;
                    nodeNew = connec(icell,2);
                    borderConnec(ielem,1) = nodeOld;
                    borderConnec(ielem,2) = nodeNew;
                end                
                
            end
            for inode = 1:nNode
                nodes = borderConnec(:,inode);
                [~,I] = sort(nodes);
                borderConnecOrdered(I,inode) = 1:length(nodes);
            end
            borderConnecSwitch(:,1) = borderConnecOrdered(:,2);
            borderConnecSwitch(:,2) = borderConnecOrdered(:,1);
        end
        
        function computeExteriorMeshesFromBoxSides(obj)
            nSides = 2;
            for iDime = 1:obj.ndim
                for iSide = 1:nSides
                    iFace = obj.computeIface(iSide,iDime);
                    [mesh,nodesInBoxFace] = obj.createBoxFaceMesh(iDime,iSide);
                    obj.boxFaceMeshes{iFace} = mesh;
                    obj.append(mesh);
                    obj.nodesInBoxFaces{iFace} = nodesInBoxFace;
                    obj.computeGlobalConnectivities(iFace);
                end
            end
            obj.nBoxFaces = numel(obj.boxFaceMeshes);                              
        end
        
        function connec = computeGlobalConnectivities(obj,iFace)
            boxFaceMesh    = obj.boxFaceMeshes{iFace};
            nodesInBoxFace = obj.nodesInBoxFaces{iFace};
            nodes = find(nodesInBoxFace);
            boxFaceConnec = boxFaceMesh.connec;
            for inode = 1:size(boxFaceConnec,2)
               connec(:,inode) = nodes(boxFaceConnec(:,inode));
            end
            obj.globalConnectivities{iFace} = connec;
        end
        
        function [m, nodesInBoxFace] = createBoxFaceMesh(obj,idime,iside)
            [boxFaceCoords,nodesInBoxFace,boxFaceCoordsRemoved] = obj.getFaceCoordinates(idime,iside);
            switch obj.ndim
                case 2
                    boxFaceConnec = obj.computeConnectivities(boxFaceCoordsRemoved);
                case 3
                    boxFaceConnec = obj.computeDelaunay(boxFaceCoordsRemoved);
            end
            s.connec = boxFaceConnec;
            s.coord  = boxFaceCoords;
            s.isInBoundary = true;
            m = Mesh().create(s);
        end
        
        function [boxFaceCoords, nodesInBoxFace,boxFaceCoordsRemoved] = getFaceCoordinates(obj,idime,iside)
            D = obj.getFaceCharacteristicDimension(idime,iside);
            nodesInBoxFace = obj.innerMeshOLD.coord(:,idime) == D;
            boxFaceCoords = obj.innerMeshOLD.coord(nodesInBoxFace,:);
            boxFaceCoordsRemoved = obj.removeExtraDimension(boxFaceCoords,idime);
            obj.storeRemovedDimensions(idime,iside,D);
        end
        
        function D = getFaceCharacteristicDimension(obj,idime,iside)
            if iside == 1
                D = min(obj.innerMeshOLD.coord(:,idime));
            elseif iside == 2
                D = max(obj.innerMeshOLD.coord(:,idime));
            else
                error('Invalid iside value. Valid values: 1 and 2.')
            end
        end
        
        function face_coord = removeExtraDimension(obj,face_coord,idime)
            dimen = 1:obj.ndim;
            face_coord = face_coord(:,dimen(dimen~=idime));
        end
        
        function storeRemovedDimensions(obj,iDime,iSide,D)
            iFace = obj.computeIface(iSide,iDime);
            obj.removedDimensions(iFace) = iDime;
            obj.removedDimensionCoord(iFace) = D;
        end
        
        function defineActiveMeshes(obj)
            obj.activeMeshesList = find([false true(1,obj.nBoxFaces)]);
            obj.nActiveMeshes     = numel(obj.activeMeshesList);
        end        
        
    end
    
    methods (Access = private, Static)
        
        function iFace = computeIface(iSide,iDime)
            nSides = 2;
            iFace = (iDime-1)*nSides + iSide;
        end        
        
        function connec = computeDelaunay(coord)
            DT = delaunayTriangulation(coord);
            connec = DT.ConnectivityList;
        end
        
        function connec = computeConnectivities(coord)
            [~,I] = sort(coord);
            connec = [I circshift(I,-1)];
            connec(end,:) = [];
        end
        
    end
    
end