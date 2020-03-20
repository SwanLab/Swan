classdef SubCellConnecComputer < handle
    
    properties (GetAccess = public, SetAccess = private)
        connec
    end
    
    properties (Access = private)
        coord
        nElem
        connecTotal
        elemCases
        
    end
    
    
    methods (Access = public)
        
        function obj = SubCellConnecComputer(cParams)
            obj.init(cParams);
            obj.compute();
        end
        
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.connecTotal = cParams.connecTotal;
            obj.nElem = cParams.nElem;
            obj.elemCases = cParams.elemCases;
            obj.coord = cParams.coord;            
        end
        
        function compute(obj)
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
                    node = localTriangleCon(icase,inode);
                    connecT(inode,elemCase) = node ;
                    allConnec(1,inode,elemCase) = connecTot(elemCase,node);
                    allConnec1(1,inode,elemCase) = connecTot(elemCase,node);
                    allConnec2(1,inode,elemCase) = connecTot(elemCase,node);
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
            
            
            
            q1 = obj.computeQuality(allConnec1,coordT);
            q2 = obj.computeQuality(allConnec2,coordT);
            
            qT(1,:) = sum(q1,1);
            qT(2,:) = sum(q2,1);
            
            [~,imax] = max(qT);
            
            subcase = imax == 1;
            allConnec(2:3,:,subcase) = allConnec1(2:3,:,subcase);
            subcase = imax == 2;
            allConnec(2:3,:,subcase) = allConnec2(2:3,:,subcase);
            obj.connec = allConnec;
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
