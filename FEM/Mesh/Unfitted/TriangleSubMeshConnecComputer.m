classdef TriangleSubMeshConnecComputer < handle
    
    properties (Access = public)
        nSubCellNodes
        nSubCellsByElem
    end
    
    properties (Access = private)
        localTriangleConnecCases
        localQuadConnecCases
        
        nodeCaseByElem
        
        imax
        cellNodes
        nSubCellsByQuad
        nSubCases
        nElemInCase  
        
        bestSubCellCaseSelectorParams
    end
    
    methods (Access = public)
        
        function obj = TriangleSubMeshConnecComputer(cParams)
            obj.init(cParams)
        end
        
        function computeBestQuadrilateralSubCell(obj,icase)
            localConnec = obj.localQuadConnecCases(:,:,:,icase);
            connecCases = obj.computeSubCasesConnec(localConnec);            
            
            s = obj.bestSubCellCaseSelectorParams;
            s.nodesSubCases = connecCases;
            bestCase = BestSubCellCaseSelector(s);
            obj.imax = bestCase.compute();            
        end
        
        function compute(obj,nodes,icase)
            obj.nElemInCase = size(nodes,1);
            obj.cellNodes   = nodes;
                        
            obj.computeBestQuadrilateralSubCell(icase);
            
            
            nodeCaseByElem = zeros(obj.nSubCellsByElem,obj.nSubCellNodes,obj.nElemInCase);
            for isubCase = 1:2
                isActive = obj.imax == isubCase;
                nodesT = obj.localTriangleConnecCases(icase,:);
                nodesQ = obj.localQuadConnecCases(:,:,isubCase,icase);
                
                for inode = 1:3                
                    nodeCaseByElem(1,inode,isActive) = nodesT(inode);
                    nodeCaseByElem(2,inode,isActive) = nodesQ(1,inode);
                    nodeCaseByElem(3,inode,isActive) = nodesQ(2,inode);
                end
            end
            
            obj.nodeCaseByElem = nodeCaseByElem;
            
        end
        
        function connec = partition(obj,nodes)
            t = obj.nodeCaseByElem;
            nElem = obj.nElemInCase;
            nnodesSubCell = size(nodes,2);
            connec = obj.initSubCellsConnec();
            for isubcell = 1:3
                for inode = 1:3
                    node = squeeze(t(isubcell,inode,:));
                    ind = sub2ind([nElem,nnodesSubCell],(1:nElem)',node);
                    connec(isubcell,inode,:) = nodes(ind);
                end
            end
            
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            s = cParams.bestSubCellCaseSelector;            
            obj.bestSubCellCaseSelectorParams = s;
            obj.nSubCellNodes   = 3;
            obj.nSubCellsByQuad = 2;
            obj.nSubCases       = 2;
            obj.nSubCellsByElem = 3;
            obj.computeTriangleLocalConnecCases();
            obj.computeLocalQuadConnecCases();
        end
        
        function computeTriangleLocalConnecCases(obj)
            nodes = [1 4 5;4 2 5; 4 3 5];
            obj.localTriangleConnecCases = nodes;
        end
        
        function computeLocalQuadConnecCases(obj)
            nodes(:,:,1,1) = [4 2 3;5 4 3];
            nodes(:,:,2,1) = [2 3 5;4 2 5];
            
            nodes(:,:,1,2) = [4 5 3;1 4 3];
            nodes(:,:,2,2) = [1 4 5;1 5 3];
            
            nodes(:,:,1,3) = [1 2 5;2 4 5];
            nodes(:,:,2,3) = [1 4 5;1 2 4];            
            obj.localQuadConnecCases = nodes;
        end
        
        function c = initSubCellsConnec(obj)
            c = zeros(obj.nSubCellsByElem,obj.nSubCellNodes,obj.nElemInCase);             
        end
        
        function connec = computeSubTriangleConnec(obj,icase)
            localConnec = obj.localTriangleConnecCases(icase,:);
            connec = zeros(obj.nSubCellNodes,obj.nElemInCase);
            for inode = 1:obj.nSubCellNodes
                localNode = localConnec(inode);
                node = obj.cellNodes(:,localNode);
                connec(inode,:) = node;
            end
            obj.connecSubTriangle = connec;            
        end
        
        function computeSubQuadConnec(obj,icase)
            localConnec = obj.localQuadConnecCases(:,:,:,icase);
            connecCases = obj.computeSubCasesConnec(localConnec);
            connec = obj.computeBestCase(connecCases);
            obj.connecSubQuad = connec;             
        end        
        
        function nodeQ = computeBestCase(obj,nodesSubCases)
           s = obj.bestSubCellCaseSelectorParams;
           s.nodesSubCases = nodesSubCases;
           bestCase = BestSubCellCaseSelector(s);
           nodeQ = bestCase.compute();
           obj.imax = bestCase.imax;
        end
        
        function connec = computeSubCasesConnec(obj,localConnecCases)
            connec = obj.initQuadConnecCases();
            for isubCase = 1:obj.nSubCases
                localConnec = localConnecCases(:,:,isubCase);
                connecCase  = obj.computeOneCaseQuadConnec(localConnec);
                connec(:,:,isubCase,:) = connecCase;
            end
        end
        
        function connec = computeOneCaseQuadConnec(obj,localConnec)
            connec = obj.initConnecQuad();
            for isubCell = 1:obj.nSubCellsByQuad
                for inode = 1:obj.nSubCellNodes
                    localNode = localConnec(isubCell,inode);
                    node = obj.cellNodes(:,localNode);
                    connec(isubCell,inode,:) = node;
                end
            end
        end
        
        function connec = initConnecQuad(obj)
            n1 = obj.nSubCellsByQuad;
            n2 = obj.nSubCellNodes;
            n3 = obj.nElemInCase;            
            connec = zeros(n1,n2,n3);
        end
        
        function connec = initQuadConnecCases(obj)
            n1 = obj.nSubCellsByQuad;
            n2 = obj.nSubCellNodes;
            n3 = obj.nSubCases;
            n4 = obj.nElemInCase;
            connec = zeros(n1,n2,n3,n4);
        end
        
    end
    
end