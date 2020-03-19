classdef CutMeshComputerProvisional < handle
    
    properties (Access = public)
        connecCutInt
        coord
    end
    
    properties (Access = private)
        backgroundConnec
        backgroundCoord
        levelSet
        
        cutEdges
        nCutEdges
        
        edgesComputer
        
        xCutPoints
        
        cutMeshPoints
        nElem
        
        isEdgeCutInElem
        cutEdgeInElem
    end
    
    methods (Access = public)
        
        function  obj = CutMeshComputerProvisional(cParams)
            obj.init(cParams);
            obj.compute()
        end
        
    end
    
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.backgroundConnec   = cParams.connec;
            obj.backgroundCoord    = cParams.coord;
            obj.levelSet = cParams.levelSet;
            obj.nElem = size(obj.backgroundConnec,1);
        end
        
        function compute(obj)
            obj.computeEdges();
            obj.computeCutEdges();
            obj.computeNcutEdges();
            obj.computeCutPoints();
            obj.computeCutMeshCoordinates();
            
            obj.computeIsEdgeCutInElem();
            obj.computeCutEdgeInElem();
            
            
            [isoNode,cutNodePerElemen] = obj.computeNodeInElementForSubTriangle(obj.backgroundCoord);
            
            isoNodeIsFull = obj.levelSet(isoNode) < 0;
            
            
            [allConnec,isSubCellInterior] = obj.createConnec(isoNodeIsFull,cutNodePerElemen);
            
            obj.connecCutInt = obj.computeCutInteriorConnec(isSubCellInterior,allConnec);
        end
        
        function computeCutMeshCoordinates(obj)
            obj.coord = [obj.backgroundCoord;obj.xCutPoints];
        end
        
        function computeEdges(obj)
            s.nodesByElem = obj.backgroundConnec;
            e = EdgesConnectivitiesComputer(s);
            e.compute();
            obj.edgesComputer = e;
        end
        
        function connecCutInterior = computeCutInteriorConnec(obj,isSubCellInterior,allConnec)
            connecTall = reshape(permute(allConnec,[1 3 2]),30,3);
            isSubCellInteriorAll = isSubCellInterior(:);
            nSubCellInerior = sum(isSubCellInteriorAll);
            connecCutInterior = zeros(nSubCellInerior,3);
            connecCutInterior(:,1) =  connecTall(isSubCellInteriorAll,1);
            connecCutInterior(:,2) =  connecTall(isSubCellInteriorAll,2);
            connecCutInterior(:,3) =  connecTall(isSubCellInteriorAll,3);
        end
        
        function [allConnec,isSubCellInterior] = createConnec(obj,isoNodeIsFull,cutNodePerElemen)
            coordT = obj.coord;
            
            connecTot = obj.computeConnecTotal(obj.backgroundCoord,obj.backgroundConnec,cutNodePerElemen);
            elemActive  = obj.isEdgeCutInElem;
            
            
            codeCase = zeros(size(elemActive,1),1);
            for iedge = 1:3
                codeCase = codeCase + elemActive(:,iedge)*2^(3 - iedge);
            end
            
            code = [5 6 3];
            nCutCells = size(codeCase,1);
            cutElemCase = zeros(size(codeCase));
            connecT = zeros(3,nCutCells);
            connecQ = zeros(4,size(codeCase,1));
            
            nSubElem = 3;
            allConnec = zeros(nSubElem,3,size(codeCase,1));
            
            localTriangleCon = [1 4 5;
                4 2 5;
                4 3 5];
            
            localQuadConnec = [4 2 3 5;
                1 4 5 3;
                1 2 4 5];
            
            
            nCases = 3;
            nNode = 2;
            nEdge = 2;
            largCase = zeros(nEdge,nNode,nCases);
            largCase(:,:,1) = [4 3; 5 2];
            
            
            largCase(:,:,2) = [4 3; 1 5];
            
            
            largCase(:,:,3) = [1 4; 5 2];
            
            
            connecCase(:,:,1,1) = [4 2 3;5 4 3];
            connecCase(:,:,1,2) = [2 3 5;4 2 5];
            
            connecCase(:,:,2,1) = [4 5 3;1 4 3];
            connecCase(:,:,2,2) = [1 4 5;1 5 3];
            
            connecCase(:,:,3,1) = [1 2 5;2 4 5];
            connecCase(:,:,3,2) = [1 4 5;1 2 4];
            
            nNodeT = 3;
            nNodeQ = 4;
            
            allConnec1 = allConnec;
            allConnec2 = allConnec;
            
            for icase = 1:nCases
                elemCases = codeCase == code(icase);
                cutElemCase(elemCases,1) = icase;
                for inode = 1:nNodeT
                    node = localTriangleCon(icase,inode);
                    connecT(inode,elemCases) = node ;
                    allConnec(1,inode,elemCases) = connecTot(elemCases,node);
                    allConnec1(1,inode,elemCases) = connecTot(elemCases,node);
                    allConnec2(1,inode,elemCases) = connecTot(elemCases,node);
                end
                
                for inode = 1:nNodeQ
                    connecQ(inode,elemCases) = localQuadConnec(icase,inode);
                end
                
                lCase = largCase(:,:,icase);
                
                nCutElem = sum(elemCases);
                nEdge = 2;
                norm = zeros(nCutElem,2);
                
                for iedge = 1:nEdge
                    nodesEdge = lCase(iedge,:);
                    x = zeros(nCutElem,2);
                    y = zeros(nCutElem,2);
                    for inode = 1:2
                        nodes = connecTot(elemCases,nodesEdge(inode));
                        x(:,inode) = coordT(nodes,1);
                        y(:,inode) = coordT(nodes,2);
                    end
                    norm(:,iedge) = sqrt((x(:,1)-x(:,2)).^2 + (y(:,1) - y(:,2)).^2);
                end
                
                [a,imax] = min(norm,[],2);
                
                subElem = find(elemCases);
                
                
                isubcase = 1;
                for inode = 1:3
                    for isubElem = 1:2
                        node = connecCase(isubElem,inode,icase,isubcase);
                        allConnec1(1+isubElem,inode,subElem) = connecTot(subElem,node);
                    end
                end
                
                
                isubcase = 2;
                for inode = 1:3
                    for isubElem = 1:2
                        node = connecCase(isubElem,inode,icase,isubcase);
                        allConnec2(1+isubElem,inode,subElem) = connecTot(subElem,node);
                    end
                end
                %end
                
                
            end
            
            
            q1 = obj.computeQuality(allConnec1,codeCase,coordT);
            q2 = obj.computeQuality(allConnec2,codeCase,coordT);
            
            qT(1,:) = sum(q1,1);
            qT(2,:) = sum(q2,1);
            
            [~,imax] = max(qT);
            
            subcase = imax == 1;
            allConnec(2:3,:,subcase) = allConnec1(2:3,:,subcase);
            subcase = imax == 2;
            allConnec(2:3,:,subcase) = allConnec2(2:3,:,subcase);
            
            isSubCellInterior = false(nSubElem,nCutCells);
            
            isSubCellInterior(1,isoNodeIsFull) = true;
            isSubCellInterior(2,isoNodeIsFull) = false;
            isSubCellInterior(3,isoNodeIsFull) = false;
            
            isSubCellInterior(1,~isoNodeIsFull) = false;
            isSubCellInterior(2,~isoNodeIsFull) = true;
            isSubCellInterior(3,~isoNodeIsFull) = true;
            
            
        end
        
        function computeCutEdges(obj)
            nodes = obj.edgesComputer.nodesInEdges;
            nodes1 = nodes(:,1);
            nodes2 = nodes(:,2);
            ls1 = obj.levelSet(nodes1);
            ls2 = obj.levelSet(nodes2);
            obj.cutEdges = xor(ls1<0,ls2<0);
        end
        
        function computeNcutEdges(obj)
            obj.nCutEdges = sum(obj.cutEdges);
        end
        
        function xCut = computeCutPoints(obj)
            nodes = obj.edgesComputer.nodesInEdges;
            nodesInCutEdges = nodes(obj.cutEdges,:);
            node1 = nodesInCutEdges(:,1);
            node2 = nodesInCutEdges(:,2);
            ls1 = obj.levelSet(node1);
            ls2 = obj.levelSet(node2);
            x1  = obj.backgroundCoord(node1,:);
            x2  = obj.backgroundCoord(node2,:);
            xCut = x1+ls1.*(x2-x1)./(ls1-ls2);
            obj.xCutPoints = xCut;
        end
        
        function isEdgeCut = computeIsEdgeCutInElem(obj)
            edgesInElem = obj.edgesComputer.edgesInElem;            
            nEdgeByElem = obj.edgesComputer.nEdgeByElem;
            isEdgeCut = false(obj.nElem,nEdgeByElem);
            for iedge = 1:nEdgeByElem
                edge = edgesInElem(:,iedge);
                isEdgeCut(:,iedge) = obj.cutEdges(edge);
            end   
            obj.isEdgeCutInElem = isEdgeCut;            
        end
        
        
        function computeCutEdgeInElem(obj)
            nCutElem = obj.nElem;
            acc = zeros(nCutElem,1);
            cutEdge = zeros(nCutElem,2);
            trow(:,1) = 1:nCutElem;
            for i = 1:3
                row = obj.isEdgeCutInElem(:,i);
                acc(row) = acc(row) + 1;
                colum = acc(row);
                A = [trow(row),colum,i*ones(size(colum))];
                tt = accumarray(A(:,1:2),A(:,3),[nCutElem 2]);
                cutEdge = cutEdge + tt;
            end
            obj.cutEdgeInElem = cutEdge;
        end
        
        function [isoNode, cutNodePerElemen] = computeNodeInElementForSubTriangle(obj,coord)
            
            eT = obj.cutEdgeInElem;
            edgesElem =  obj.edgesComputer.edgesInElem;
            edgesNodes = obj.edgesComputer.nodesInEdges;            
            edgesInCutElemens = edgesElem(:,:);
            
            [rows,cols] = find(obj.isEdgeCutInElem);
            nodesWithCtuEdge = zeros(size(eT,1),4);
            
            nCutElem = obj.nElem;
            nIn = size(coord,1);
            
            
            cutNodes = zeros(size(edgesNodes,1),1);
            cutNodes(obj.cutEdges) = (nIn+1):(nIn+obj.nCutEdges);
            
            nCutEdges = 2;
            cutNodePerElemen = zeros(nCutElem,nCutEdges);
            
            for it = 1:nCutEdges
                
                t = sub2ind([10 3],[1:10]',eT(:,it));
                edge = edgesInCutElemens(t);
                
                cutNodePerElemen(:,it) = cutNodes(edge);
                
                sEdge = edgesNodes(edge,:);
                for inode = 1:2
                    node = sEdge(:,inode);
                    ct = 2*(it -1) + inode;
                    nodesWithCutEdge(:,ct) =  node;
                    
                end
            end
            
            [isoNode,b] = mode(nodesWithCutEdge,2);
            
            
        end
        
        %
        function connecT = computeConnecTotal(obj,coord,connec,cutNodePerElemen)
            edgesElem = obj.edgesComputer.edgesInElem;                        
            connecCut = connec(:,:);
            
            edgesCutElem = edgesElem(:,:);
            
            nIn = size(coord,1);
            nCutNode = 2;
            nnode = size(connec,2);
            nNodeConnec = nnode + nCutNode;
            nCutElem = obj.nElem;
            connecT = zeros(size(connecCut,1),nNodeConnec);
            nodesCut(:,1) = (nIn+1):(nIn + nCutElem);
            
            connecT(:,1) = connecCut(:,1);
            connecT(:,2) = connecCut(:,2);
            connecT(:,3) = connecCut(:,3);
            
            connecT(:,4) = cutNodePerElemen(:,1);
            connecT(:,5) = cutNodePerElemen(:,2);
            
            
        end
        
    end
        
        methods (Access = private, Static)
        
        
        function q = computeQuality(allConnec,codeCase,coordT)
            A = zeros(2,size(codeCase,1));
            q = zeros(2,size(codeCase,1));
            for isubElem = 1:2
                
                nodes = squeeze(allConnec(1+isubElem,:,:))';
                
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