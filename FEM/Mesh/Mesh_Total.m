classdef Mesh_Total < Mesh_Composite 
    
    properties (GetAccess = public, SetAccess = private)
        nBoxFaces
        
        innerMeshOLD
        boxFaceMeshes
        nodesInBoxFaces
        

        
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
            obj.type = obj.innerMeshOLD.type;
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
            obj.innerMeshOLD = Mesh(s);
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
                s.coord  = obj.coord(nodes,:);
                s.connec = obj.computeConnectivitiesFromData();
                %s.globalConnec = obj.borderElements(:,2:end);
                s.nodesInBoxFaces = false(size(obj.coord,1),1);
                s.nodesInBoxFaces(nodes,1) = true;
                m = BoundaryMesh(s);
                obj.boxFaceMeshes{imesh} = m;
                obj.append(m);
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
               s.backgroundMesh = obj.innerMeshOLD;
               s.dimensions = 1:s.backgroundMesh.ndim;
               bC = BoundaryMeshCreatorFromRectangularBox(s);
               bMeshes = bC.create();
               obj.nBoxFaces = numel(bMeshes);
               for iM = 1:obj.nBoxFaces
                  m = bMeshes{iM};
                  obj.boxFaceMeshes{iM} = m;
                  obj.append(m);                   
               end
        end
        
        function defineActiveMeshes(obj)
            obj.activeMeshesList = find([false true(1,obj.nBoxFaces)]);
            obj.nActiveMeshes     = numel(obj.activeMeshesList);
        end        
        
    end
    
    
end