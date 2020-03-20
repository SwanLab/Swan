classdef CutMeshComputerProvisional < handle
    
    properties (Access = public)
        connec
        coord
    end
    
    properties (Access = private)
        backgroundConnec
        backgroundCoord
        levelSet
        
        edgesComputer
        cutEdgesComputer        
        
        xCutPoints
        
        nElem
        
   
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

            obj.computeCutPoints();
            obj.computeCutMeshCoordinates();                            
          
            obj.computeCutInteriorConnec();
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
        
        function computeCutEdges(obj)
            s.levelSet = obj.levelSet;
            s.edgesComputer = obj.edgesComputer;
            c = CutEdgesComputer(s);
            c.compute();
            obj.cutEdgesComputer = c;    
        end
        
        function computeCutInteriorConnec(obj)
            isSubCellInterior = obj.computeIsSubCellInterior();  
            [allConnec] = obj.createConnec();
            
            connecTall = reshape(permute(allConnec,[1 3 2]),30,3);
            isSubCellInterior = isSubCellInterior(:);
            nSubCellInerior = sum(isSubCellInterior);
            connecT = zeros(nSubCellInerior,3);
            connecT(:,1) =  connecTall(isSubCellInterior,1);
            connecT(:,2) =  connecTall(isSubCellInterior,2);
            connecT(:,3) =  connecTall(isSubCellInterior,3);
            
            obj.connec = connecT;
        end
        
        function isSubCellInterior = computeIsSubCellInterior(obj)
            isoNode = obj.computeIsoNode();                        
            isoNodeIsFull = obj.levelSet(isoNode) < 0;            
            nSubElem = 3;
            isSubCellInterior = false(nSubElem,obj.nElem);
            
            isSubCellInterior(1,isoNodeIsFull) = true;
            isSubCellInterior(2,isoNodeIsFull) = false;
            isSubCellInterior(3,isoNodeIsFull) = false;
            
            isSubCellInterior(1,~isoNodeIsFull) = false;
            isSubCellInterior(2,~isoNodeIsFull) = true;
            isSubCellInterior(3,~isoNodeIsFull) = true;            
        end        
        
        function [allConnec] = createConnec(obj)
            s.connecTotal = obj.computeConnecTotal();
            s.nElem = obj.nElem;
            s.elemCases = obj.cutEdgesComputer.elemCases;
            s.coord = obj.coord;
            subCell = SubCellConnecComputer(s);
            allConnec = subCell.connec;
        end
        
        function xCut = computeCutPoints(obj)
            nodes = obj.edgesComputer.nodesInEdges;
            cutEdges = obj.cutEdgesComputer.cutEdges;
            nodesInCutEdges = nodes(cutEdges,:);
            node1 = nodesInCutEdges(:,1);
            node2 = nodesInCutEdges(:,2);
            ls1 = obj.levelSet(node1);
            ls2 = obj.levelSet(node2);
            x1  = obj.backgroundCoord(node1,:);
            x2  = obj.backgroundCoord(node2,:);
            xCut = x1+ls1.*(x2-x1)./(ls1-ls2);
            obj.xCutPoints = xCut;
        end
      
        
        
            
        function isoNode = computeIsoNode(obj)
            nEdges = obj.cutEdgesComputer.nCutEdgeByElem;            
            cutEdgeInElem = obj.cutEdgesComputer.cutEdgeInElem;            
            edgesNodes = obj.edgesComputer.nodesInEdges;                  
            for iedge = 1:nEdges
                edge = cutEdgeInElem(:,iedge);
                sEdge = edgesNodes(edge,:);
                for inode = 1:2
                    node = sEdge(:,inode);
                    ct = 2*(iedge -1) + inode;
                    nodesWithCutEdge(:,ct) =  node;                   
                end
            end                                   
            [isoNode,~] = mode(nodesWithCutEdge,2);
        end
        
        %
        function connecT = computeConnecTotal(obj)
            
         
            cutNodePerElemen = obj.computeCutNodePerElem();            
            
            nCutNode = obj.cutEdgesComputer.nCutEdgeByElem;
            nnode    = size(obj.backgroundConnec,2);
           
            nNodeConnec = nnode + nCutNode;

            bConnec = obj.backgroundConnec;            
            connecT = zeros(obj.nElem,nNodeConnec);
            
            connecT(:,1) = bConnec(:,1);
            connecT(:,2) = bConnec(:,2);
            connecT(:,3) = bConnec(:,3);
            
            connecT(:,4) = cutNodePerElemen(:,1);
            connecT(:,5) = cutNodePerElemen(:,2);
            
            
        end
        
        function cutNodePerElemen = computeCutNodePerElem(obj)
            firstCutEdgePerElem = obj.cutEdgesComputer.firstCutEdge;
            finalNode = size(obj.backgroundCoord,1);                                    
            cutNodePerElemen = firstCutEdgePerElem + finalNode;
        end        
        
   
        
    
    end
    
    
    
end