classdef BoundaryMeshCreatorFromData < BoundaryMeshCreator
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
       borderNodes
       borderElements
       backgroundMesh
    end
    
    methods (Access = public)
        
        function obj = BoundaryMeshCreatorFromData(cParams)
            obj.init(cParams)
        end
        
        function m = create(obj)
            nodes    = obj.borderNodes;
            s.coord  = obj.backgroundMesh.coord(nodes,:);
            kConnec = boundary(s.coord);
            s.connec = [kConnec(1:end-1),kConnec(2:end)];
            %s.connec = obj.computeConnectivitiesFromData(obj.borderElements(:,2:end));
           
            s.nodesInBoxFaces = false(size(obj.backgroundMesh.coord,1),1);
            s.nodesInBoxFaces(nodes,1) = true;
            s.isRectangularBox = false;
            s.dimension = 1;
            s.kFace      = obj.backgroundMesh.kFace;
            m{1} = BoundaryMesh(s);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.borderNodes    = cParams.borderNodes;
            obj.borderElements = cParams.borderElements;
            obj.backgroundMesh = cParams.backgroundMesh;
        end
        
        function borderConnecSwitch = computeConnectivitiesFromData(obj,connec)
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
        
    end
    
end