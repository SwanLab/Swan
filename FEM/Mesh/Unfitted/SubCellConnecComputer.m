classdef SubCellConnecComputer < handle
    
    properties (GetAccess = public, SetAccess = private)
        connec
    end
    
    properties (Access = private)
        coord
        nElem
        connecTotal
        connecTall
        elemCases
        levelSet
        isoNode
        vertexNodesInElem
        cutNodesInElem
        allNodesInElem
        
        
        isSubCellInterior
    end
    
    methods (Access = public)
        
        function obj = SubCellConnecComputer(cParams)
            obj.init(cParams);
            obj.compute();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.cutNodesInElem   = cParams.cutNodePerElemen;
            obj.vertexNodesInElem = cParams.backgroundConnec;
            obj.nElem            = cParams.nElem;
            obj.elemCases        = cParams.elemCases;
            obj.coord            = cParams.coord;
            obj.levelSet         = cParams.levelSet;
        end
        
        function compute(obj)
            obj.computeAllNodesInElem();
            obj.computeAllConnec();
            obj.computeIsoNode();
            obj.computeIsSubCellInterior();
            obj.computeConnec();
        end
        
        function computeAllNodesInElem(obj)
            vertexNodes = obj.vertexNodesInElem;
            cutNodes    = obj.cutNodesInElem;
            allNodes    = [vertexNodes,cutNodes];
            obj.allNodesInElem = allNodes;
        end
        
        function computeAllConnec(obj)
            subTriangleConnec = obj.computeSubTriangleConnec();
            allConnecQuad     = obj.computeSubTriangleConnecFromQuad();
            
            nSubCell = 3;
            nnode    = 3;
            allConnec = zeros(nSubCell,nnode,obj.nElem);
            allConnec(1,:,:) = subTriangleConnec;
            allConnec(2:3,:,:) = allConnecQuad;
            
            allConnec = permute(allConnec,[1 3 2]);
            
            connecTall = reshape(allConnec,30,3);
            
            obj.connecTall = connecTall;
        end
        
        function allConnec = computeSubTriangleConnec(obj)
            nCases = 3;
            nNodeT = 3;
            nCutCells = obj.nElem;
            connecTot = obj.allNodesInElem;
            
            localTriangleCon = [1 4 5;
                4 2 5;
                4 3 5];
            
            connecT = zeros(3,nCutCells);
            allConnec = zeros(3,obj.nElem);
            
            for icase = 1:nCases
                elemCase = obj.elemCases(:,icase);
                for inode = 1:nNodeT
                    localNode = localTriangleCon(icase,inode);
                    connecT(inode,elemCase) = localNode ;
                    node = connecTot(elemCase,localNode);
                    allConnec(inode,elemCase) = node;
                end
                
            end
        end
        
        function allConnecQuad = computeSubTriangleConnecFromQuad(obj)            
            allConnec = obj.computeConnecSubCases();
            imax = obj.computeBetterSubCaseOption(allConnec);
            allConnecQuad = obj.computeAllConnecQuad(allConnec,imax);
        end        
        
        function allConnec = computeConnecSubCases(obj)
            connecCase(:,:,1,1) = [4 2 3;5 4 3];
            connecCase(:,:,1,2) = [2 3 5;4 2 5];
            
            connecCase(:,:,2,1) = [4 5 3;1 4 3];
            connecCase(:,:,2,2) = [1 4 5;1 5 3];
            
            connecCase(:,:,3,1) = [1 2 5;2 4 5];
            connecCase(:,:,3,2) = [1 4 5;1 2 4];
            
            
            
            allConnec = zeros(2,3,obj.nElem,2);
            
            for isubCase = 1:2
                connecT = connecCase(:,:,:,isubCase);
                connecSubCase = obj.computeConnecSubCase(connecT);
                allConnec(:,:,:,isubCase) = connecSubCase;
            end
            
            
        end
        
        
        function allConnec = computeConnecSubCase(obj,connec)
            allConnec = zeros(2,3,obj.nElem);
            connecTot = obj.allNodesInElem;
            
            nCases = 3;
            
            for icase = 1:nCases
                elemCase = obj.elemCases(:,icase);
                connec2 = connecTot(elemCase,:);
                for inode = 1:3
                    for isubElem = 1:2
                        localNode = connec(isubElem,inode,icase);
                        node = connec2(:,localNode);
                        allConnec(isubElem,inode,elemCase) = node;
                    end
                end
                
            end
        end
        
        function imax = computeBetterSubCaseOption(obj,allConnec)
            qT = zeros(2,obj.nElem);
            for isubCase = 1:2
                q = obj.computeQuality(allConnec(:,:,:,isubCase));
                qT(isubCase,:) = sum(q,1);
            end
            [~,imax] = max(qT);            
        end

        function allConnecQuad = computeAllConnecQuad(obj,allConnec,imax)
            allConnecQuad = zeros(2,3,obj.nElem);
            for isubCase = 1:2
                subcase = imax == isubCase;
                allConnecQuad(:,:,subcase) = allConnec(:,:,subcase,isubCase);
            end            
        end
        
        function q = computeQuality(obj,allConnec)
            coordT = obj.coord;
            A = zeros(2,obj.nElem);
            q = zeros(2,obj.nElem);
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
        
        function computeIsoNode(obj)
            
            isoNodeT = zeros(obj.nElem,1);
            
            nCases = 3;
            
            isoNodeCase = [1 2 3];
            
            for icase = 1:nCases
                elemCase = obj.elemCases(:,icase);
                isoNodeT(elemCase,1) = isoNodeCase(icase);
            end
            
            t = sub2ind([10 3],[1:10]',isoNodeT);
            
            obj.isoNode = obj.allNodesInElem(t);
        end
        
        function isSubCellInterior = computeIsSubCellInterior(obj)
            
            isoNodeIsFull = obj.levelSet(obj.isoNode) < 0;
            
            isTriangleInterior = isoNodeIsFull;
            
            
            nSubElem = 3;
            isSubCellInterior = false(nSubElem,obj.nElem);
            
            isSubCellInterior(1,isTriangleInterior) = true;
            isSubCellInterior(2,isTriangleInterior) = false;
            isSubCellInterior(3,isTriangleInterior) = false;
            
            isSubCellInterior(1,~isTriangleInterior) = false;
            isSubCellInterior(2,~isTriangleInterior) = true;
            isSubCellInterior(3,~isTriangleInterior) = true;
            
            obj.isSubCellInterior = isSubCellInterior;
        end
        
        
        function computeConnec(obj)
            isInterior = obj.isSubCellInterior(:);
            nSubCellInerior = sum(isInterior);
            connecT = zeros(nSubCellInerior,3);
            connecT(:,1) =  obj.connecTall(isInterior,1);
            connecT(:,2) =  obj.connecTall(isInterior,2);
            connecT(:,3) =  obj.connecTall(isInterior,3);
            obj.connec = connecT;
        end
        
        
        
    end
    
end
