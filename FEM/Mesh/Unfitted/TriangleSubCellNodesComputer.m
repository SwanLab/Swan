classdef TriangleSubCellNodesComputer < handle
    
    properties (Access = public)
        allSubCellsConnec
        isoNode
    end
    
    properties (Access = private)
        localQuadNodeCases
        localTriangleNodeCases
        
        nSubCases
        nSubCellsByQuad
        
        nodesInElem
        nSubCellNodes
        nSubCellsByElem
        nElem
        nCases
        nSubCells        
        
        allNodesInElem
        
        elemCases
        
        coord
    end
    
    methods (Access = public)
        
        % SUBMESHERRRRR
        function obj = TriangleSubCellNodesComputer(cParams)
            obj.init(cParams)
        end
        
        function compute(obj)
            
            isoNodeCase = [1 2 3];
            
            nodesInSubCells = zeros(obj.nSubCellsByElem,obj.nSubCellNodes,obj.nElem);
            
            isoNodeT = zeros(obj.nElem,1);
            
            for icase = 1:obj.nCases
                elemCase     = obj.elemCases(:,icase);
                nElemInCase  = sum(elemCase);
                
                nodes        = obj.nodesInElem(elemCase,:);
                
                localNodes   = obj.localTriangleNodeCases(icase,:);
                nodesT = obj.computeNodesSubCasesTriangle(localNodes,nodes,nElemInCase);
                nodesInSubCells(1,:,elemCase) = nodesT;
                
                
                localNodes = obj.localQuadNodeCases(:,:,:,icase);
                nodeQuad = obj.computeAllConnecQuad(localNodes,nodes,nElemInCase);
                nodesInSubCells(2:3,:,elemCase) = nodeQuad;
                
                isoNodeT(elemCase,1) = isoNodeCase(icase);
                
            end
            obj.allSubCellsConnec = obj.computeAllSubCells(nodesInSubCells);
            
            t = sub2ind([obj.nElem obj.nSubCellsByElem],(1:obj.nElem)',isoNodeT);
            obj.isoNode = obj.nodesInElem(t);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.nSubCases       = 2;
            obj.nSubCellsByQuad = 2;
            
            obj.coord = cParams.coord;
        
            obj.nodesInElem = cParams.nodesInElem;
            obj.nSubCellNodes = cParams.nSubCellNodes;
            obj.nSubCellsByElem = cParams.nSubCellsByElem;
            obj.nElem = cParams.nElem;
            obj.nCases = cParams.nCases;                
            obj.elemCases = cParams.elemCases;   
            obj.nSubCells = cParams.nSubCells;
            
            obj.computeTriangleLocalNodeCases();
            obj.computeLocalNodesCases();
            
        end
        
        function computeTriangleLocalNodeCases(obj)
            nodes = [1 4 5;4 2 5; 4 3 5];
            obj.localTriangleNodeCases = nodes;
        end
        
        function computeLocalNodesCases(obj)
            nodes(:,:,1,1) = [4 2 3;5 4 3];
            nodes(:,:,2,1) = [2 3 5;4 2 5];
            
            nodes(:,:,1,2) = [4 5 3;1 4 3];
            nodes(:,:,2,2) = [1 4 5;1 5 3];
            
            nodes(:,:,1,3) = [1 2 5;2 4 5];
            nodes(:,:,2,3) = [1 4 5;1 2 4];
            
            obj.localQuadNodeCases = nodes;
        end
        
        function nodes = computeAllSubCells(obj,nodes)
            nodes = permute(nodes,[1 3 2]);
            nodes = reshape(nodes,obj.nSubCells,obj.nSubCellNodes);
        end
        
        
        function nodesSubCases = computeNodesSubCasesTriangle(obj,localNodes,nodes,nElemInCase)
            nodesSubCases = zeros(obj.nSubCellNodes,nElemInCase);
            for inode = 1:obj.nSubCellNodes
                localNode = localNodes(inode);
                node = nodes(:,localNode);
                nodesSubCases(inode,:) = node;
            end
        end
        
        function nodeQuad = computeAllConnecQuad(obj,localNodesCase,nodes,nElemInCase)
            nodesSubCases = obj.computeNodesSubCases(nodes,localNodesCase,nElemInCase);
            nodeQuad = zeros(obj.nSubCellsByQuad,obj.nSubCellNodes,nElemInCase);
            imax = obj.computeBetterSubCaseOption(nodesSubCases,nElemInCase);
            for isubCase = 1:obj.nSubCases
                isSubCaseActive = imax == isubCase;
                nodesSubCase = squeeze(nodesSubCases(:,:,isubCase,isSubCaseActive));
                nodeQuad(:,:,isSubCaseActive) = nodesSubCase;
            end
        end
        
        function nodesSubCases = computeNodesSubCases(obj,nodes,localNodesCase,nElemInCase)
            nodesSubCases = zeros(obj.nSubCellsByQuad,obj.nSubCellNodes,obj.nSubCases,nElemInCase);
            for isubCase = 1:obj.nSubCases
                localNodes = localNodesCase(:,:,isubCase);
                nodesT = zeros(obj.nSubCellsByQuad,obj.nSubCellNodes,nElemInCase);
                for isubCell = 1:obj.nSubCellsByQuad
                    for inode = 1:obj.nSubCellNodes
                        localNode = localNodes(isubCell,inode);
                        node = nodes(:,localNode);
                        nodesT(isubCell,inode,:) = node;
                    end
                end
                nodesSubCases(:,:,isubCase,:) = nodesT;
            end
        end
        
        function imax = computeBetterSubCaseOption(obj,nodesSubCases,nElemInCase)
            qT = zeros(obj.nSubCases,nElemInCase);
            for isubCase = 1:obj.nSubCases
                nodesT = nodesSubCases(:,:,isubCase,:);
                q = obj.computeQuality(nodesT,nElemInCase);
                qT(isubCase,:) = sum(q,1);
            end
            [~,imax] = max(qT);
        end
        
        function q = computeQuality(obj,allConnec,nElemInCase)
            coordT = obj.coord;
            A = zeros(2,nElemInCase);
            q = zeros(2,nElemInCase);
            for isubElem = 1:2
                
                nodes = squeeze(allConnec(isubElem,:,:))';
                
                xA = coordT(nodes(:,1),1);
                yA = coordT(nodes(:,1),2);
                
                xB = coordT(nodes(:,2),1);
                yB = coordT(nodes(:,2),2);
                
                xC = coordT(nodes(:,3),1);
                yC = coordT(nodes(:,3),2);
                
                A(isubElem,:) =1/2*abs((xA -xC).*(yB-yA)-(xA-xB).*(yC-yA));
                Lab = (xA - xB).^2 + (yA - yB).^2;
                Lcb = (xC - xB).^2 + (yC - yB).^2;
                Lac = (xA - xC).^2 + (yA - yC).^2;
                L = Lab + Lcb + Lac;
                q(isubElem,:) = 4*sqrt(3)*A(isubElem,:)./L';
            end
            
            
            
        end
        
    end
    
    
    
end