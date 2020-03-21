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
        backgroundConnec
        nCutEdgeByElem
        cutNodePerElemen
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
            obj.nCutEdgeByElem = cParams.nCutEdgeByElem;
            obj.cutNodePerElemen = cParams.cutNodePerElemen;
            obj.backgroundConnec = cParams.backgroundConnec;
            obj.nElem = cParams.nElem;
            obj.elemCases = cParams.elemCases;
            obj.coord = cParams.coord;
            obj.levelSet = cParams.levelSet;
        end
        
        function compute(obj)
            obj.computeConnecTotal();
            obj.computeAllConnec();
            obj.computeIsoNode();
            obj.computeIsSubCellInterior();  
            obj.computeConnec();
        end
        
        function connecT = computeConnecTotal(obj)
         
            
            nCutNode = obj.nCutEdgeByElem;
            nnode    = size(obj.backgroundConnec,2);
           
            nNodeConnec = nnode + nCutNode;

            bConnec = obj.backgroundConnec;            
            connecT = zeros(obj.nElem,nNodeConnec);
            
            connecT(:,1) = bConnec(:,1);
            connecT(:,2) = bConnec(:,2);
            connecT(:,3) = bConnec(:,3);
            
            connecT(:,4) = obj.cutNodePerElemen(:,1);
            connecT(:,5) = obj.cutNodePerElemen(:,2);
            
            obj.connecTotal = connecT;
            
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
                
        
        
        function computeAllConnec(obj)
            coordT = obj.coord;
            
            connecTot = obj.connecTotal;
            
            
            nCutCells = obj.nElem;
            connecT = zeros(3,nCutCells);
            connecQ = zeros(4,obj.nElem);
            
            nSubElem = 3;
            allConnec = zeros(nSubElem,3,obj.nElem);
            
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
                elemCase = obj.elemCases(:,icase);
                % cutElemCase(elemCases,1) = icase;
                
                for inode = 1:nNodeT
                    localNode = localTriangleCon(icase,inode);
                    connecT(inode,elemCase) = localNode ;
                    node = connecTot(elemCase,localNode);
                    allConnec(1,inode,elemCase) = node;
                    allConnec1(1,inode,elemCase) = node;
                    allConnec2(1,inode,elemCase) = node;
                end
                
                for inode = 1:nNodeQ
                    connecQ(inode,elemCase) = localQuadConnec(icase,inode);
                end
                
                lCase = largCase(:,:,icase);
                
                nCutElem = sum(elemCase);
                nEdge = 2;
                norm = zeros(nCutElem,2);
                
                for iedge = 1:nEdge
                    nodesEdge = lCase(iedge,:);
                    x = zeros(nCutElem,2);
                    y = zeros(nCutElem,2);
                    for inode = 1:2
                        nodes = connecTot(elemCase,nodesEdge(inode));
                        x(:,inode) = coordT(nodes,1);
                        y(:,inode) = coordT(nodes,2);
                    end
                    norm(:,iedge) = sqrt((x(:,1)-x(:,2)).^2 + (y(:,1) - y(:,2)).^2);
                end
                
                [a,imax] = min(norm,[],2);
                
                subElem = find(elemCase);
                
                
                isubcase = 1;
                for inode = 1:3
                    for isubElem = 1:2
                        localNode = connecCase(isubElem,inode,icase,isubcase);
                        allConnec1(1+isubElem,inode,subElem) = connecTot(subElem,localNode);
                    end
                end
                
                
                isubcase = 2;
                for inode = 1:3
                    for isubElem = 1:2
                        localNode = connecCase(isubElem,inode,icase,isubcase);
                        allConnec2(1+isubElem,inode,subElem) = connecTot(subElem,localNode);
                    end
                end
                %end
                
            end
            
            q1 = obj.computeQuality(allConnec1,coordT);
            q2 = obj.computeQuality(allConnec2,coordT);
            
            qT(1,:) = sum(q1,1);
            qT(2,:) = sum(q2,1);
            
            [~,imax] = max(qT);
            
            subcase = imax == 1;
            allConnec(2:3,:,subcase) = allConnec1(2:3,:,subcase);
            subcase = imax == 2;
            allConnec(2:3,:,subcase) = allConnec2(2:3,:,subcase);
            
            connecTall = reshape(permute(allConnec,[1 3 2]),30,3);

            obj.connecTall = connecTall;
            
            
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
            connecTot = obj.connecTotal;
            
            obj.isoNode = connecTot(t);
        end
        
        function q = computeQuality(obj,allConnec,coordT)
            A = zeros(2,obj.nElem);
            q = zeros(2,obj.nElem);
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
