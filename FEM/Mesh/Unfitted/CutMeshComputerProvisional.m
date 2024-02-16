classdef CutMeshComputerProvisional < CutMesh
  
    properties (Access = private)
        connec
        coord
        
        cutEdgesComputer
        cutPointsInElemComputer
        
        cutEdgesParams
        cutCoordParams
        interiorSubCellsParams
        cutEdgesComputerParams
        
        cutCoordComputer
        
        subCellCases
        
        isSubCellInterior
        
        connecB
    end
        
    methods (Access = public)
        
        function  obj = CutMeshComputerProvisional(cParams)
            obj.init(cParams);
            obj.computeAllParams();
        end
        
        function compute(obj)
            obj.computeSubCellCases();
            obj.computeCutEdges();
            obj.computeCutCoordinateComputer();
            obj.coord = obj.cutCoordComputer.allValues;
            obj.computeConnecNew();
            obj.computeMesh();
            obj.computeBoundaryMesh();
            obj.computeInnerCutMesh();
            obj.computeBoundaryCutMesh();
        end
        
    end
    
    methods (Access = private)
                
        function computeAllParams(obj)
            obj.backgroundMesh.computeEdges();
            e = obj.backgroundMesh.edges;
            obj.cutEdgesParams.nodesInEdges = e.nodesInEdges;
            obj.cutEdgesParams.levelSet     = obj.levelSet;
           
            obj.cutCoordParams.fValues = obj.backgroundMesh.coord;
            obj.cutCoordParams.nodesInEdges = e.nodesInEdges;
            
            cEparams = obj.cutEdgesComputerParams;
            
            cEparams.allNodesinElemParams.finalNodeNumber = size(obj.backgroundMesh.coord,1);
            cEparams.allNodesinElemParams.connec = obj.backgroundMesh.connec;
            cEparams.allNodesInElemCoordParams.localNodeByEdgeByElem = e.localNodeByEdgeByElem;
            cEparams.edgesInElem = e.edgesInElem;
            cEparams.nEdgeByElem = e.nEdgeByElem;
            
            obj.cutEdgesComputerParams = cEparams;
            

        end
        
        function computeCutEdges(obj)
            s = obj.cutEdgesParams;
            c = CutEdgesComputer(s);
            c.compute();
            obj.cutEdgesComputer = c;
        end
        
        function isEdgeCutInElem = computeIsEdgeCutInElem(obj)
            s.edgesInElem = obj.backgroundMesh.edges.edgesInElem;
            s.isEdgeCut   = obj.cutEdgesComputer.isEdgeCut;
            isEdgeCut = EdgeCutInElemComputer(s);
            isEdgeCutInElem = isEdgeCut.compute();
        end
        
        function computeSubCellCases(obj)
            s.connec = obj.backgroundMesh.connec;
            s.levelSet = obj.levelSet;
            obj.subCellCases = SubCellsCasesComputer(s);
            obj.subCellCases.compute();
        end
     
        function computeCutCoordinateComputer(obj)
            s = obj.cutCoordParams;
            s.xCutEdgePoint    = obj.cutEdgesComputer.xCutEdgePoint;
            s.isEdgeCut        = obj.cutEdgesComputer.isEdgeCut;
            cComputer = CutFunctionValuesComputer(s);
            cComputer.compute();
            obj.cutCoordComputer = cComputer;
        end
        
        function computeCutPointsInElemComputer(obj)
            isEdgeCutInElem =  obj.computeIsEdgeCutInElem();
            
            s = obj.cutEdgesComputerParams;
            s.isEdgeCut = obj.cutEdgesComputer.isEdgeCut;
            s.allNodesInElemCoordParams.xCutEdgePoint = obj.cutEdgesComputer.xCutEdgePoint;
            s.isEdgeCutInElem = isEdgeCutInElem;
            
            sA.isEdgeCutInElem = obj.computeIsEdgeCutInElem();
            s.all2Cut = AllEdges2CutEdgesComputer(sA);
                    
            
            c = CutPointsInElemComputer(s);
            c.compute();
            obj.cutPointsInElemComputer = c;
        end
         
        function computeConnec(obj)
            c = obj.cutPointsInElemComputer;
            sS.bestSubCellCaseSelector.coord = obj.coord;
            sA.subMeshConnecParams = sS;
            sA.xAllNodesInElem = c.xAllNodesInElem;
            sA.allNodesInElem  = c.allNodesInElem;
            sA.subCellCases    = obj.subCellCases;

            s.allSubCellsConnecParams = sA;
            
            s.cutElems = obj.cutCells;
            
            switch size(obj.subCellCases,2)
                case 3
                    nSubCellsByElem   = 3;
                case 4
                    switch mode(size(c.allNodesInElem,2))
                        case 7
                            s.nSubCellsByElem = 4;
                        case 8
                            s.nSubCellsByElem = 6;
                    end
                    
                    nSubCellsByElem   = 4;
            end
            
            s.nSubCellsByElem = nSubCellsByElem;
            s.isSubCellInterior = obj.isSubCellInterior;
            subCell = InteriorSubCellsConnecComputer(s);
            obj.connec                = subCell.connec;
            obj.xCoordsIso            = subCell.xCoordsIso;
            obj.cellContainingSubcell = subCell.cellContainingSubcell;
        end

        function cE = computeCutPointInElem(obj,isEdgeCutInElemCase,t)
                e = obj.backgroundMesh.edges;
            
                s.isEdgeCutInElem = isEdgeCutInElemCase;
                all2Cut = AllEdges2CutEdgesComputer(s);            
                cEp.all2Cut = all2Cut;
                cEp.allNodesinElemParams.finalNodeNumber = size(obj.backgroundMesh.coord,1);
                cEp.allNodesinElemParams.connec = obj.backgroundMesh.connec(t,:);
                cEp.allNodesInElemCoordParams.localNodeByEdgeByElem = e.localNodeByEdgeByElem(t,:,:);
                cEp.edgesInElem = e.edgesInElem(t,:);
                cEp.nEdgeByElem = e.nEdgeByElem;
                cEp.isEdgeCut = obj.cutEdgesComputer.isEdgeCut;
                cEp.allNodesInElemCoordParams.xCutEdgePoint = obj.cutEdgesComputer.xCutEdgePoint;
                cEp.isEdgeCutInElem = isEdgeCutInElemCase;
                cE = CutPointsInElemComputer(cEp);
                cE.compute();
        end
        
        function computeConnecNew(obj)
            isEdgeCutInElem =  obj.computeIsEdgeCutInElem();
            
            nEdgesCutCase   = [2 3 4];
            nSubCellsByElem = [3 4 6];
            
            subCell = cell(length(nEdgesCutCase),1);
            cN = cell(length(nEdgesCutCase),1);
            tP = cell(length(nEdgesCutCase),1);
            
            nCutEdges = sum(isEdgeCutInElem,1);
            for icases = 1:length(nEdgesCutCase)
                t = nCutEdges == nEdgesCutCase(icases);
                isEdgeCutInElemCase = isEdgeCutInElem(:,t);
                         
                cE = obj.computeCutPointInElem(isEdgeCutInElemCase,t);

                tP{icases} = t;
                
                if sum(t) ~= 0
                    
                    cN{icases} = cE;
                    
                    caseInfo = obj.subCellCases.caseInfo{icases};
                    
                    
                    nodes = obj.backgroundMesh.connec;
                    cutCells(:,1) = 1:size(nodes,1);
                    
                    
                    
                    sS.bestSubCellCaseSelector.coord = obj.coord;
                    sA.subMeshConnecParams           = sS;
                    sA.xAllNodesInElem               = cE.xAllNodesInElem;
                    sA.allNodesInElem                = cE.allNodesInElem;
                    sA.subCellCases                  = caseInfo.subCellCases(t,:);
                    
                    sI.allSubCellsConnecParams = sA;
                    sI.isSubCellInterior = caseInfo.isSubCellsInterior(:,t);
                    sI.cutElems = obj.cutCells(t);
                    
                    sI.nSubCellsByElem = nSubCellsByElem(icases);
                    
                    subCell{icases} = InteriorSubCellsConnecComputer(sI);
                end
                
            end
            
            connecT = [];
            xCoordsIso = [];
            cellC = [];
            
            connecBound = [];
            xCoordsIsoBoundary = [];
            cellCB = [];
            
            for icase = 1:3
                subC = subCell{icase};
                scN  = cN{icase};
                if ~isempty(subC)
                    connecT = cat(1,connecT,subC.connec);
                    xCoordsIso = cat(3,xCoordsIso,subC.xCoordsIso);
                    cellC = cat(1,cellC,subC.cellContainingSubcell);
                    
                    
                    
                end
                
                if ~isempty(scN)
                    conn = scN.cutNodesInElem;
                    xB = scN.xCutInElem;
                    tCC = tP{icase};
                    localIncidenceMatrix = scN.nodesInCutEdges;
                    
                    cellCBA = obj.cutCells(tCC);
                    if size(conn,2) == 4
                        
                        s.localIncidenceMatrix = localIncidenceMatrix;
                        s.coord        = obj.coord;
                        s.connec       = conn;
                        s.xLocalVertex = xB;
                        qS = QuadrilaterSubTriangulator(s);
                        [connR,xCoordsIsoBoundaryR] = qS.compute();
                        
                        cellCBR = [cellCBA;cellCBA];
                        
                    else
                        connR = conn;
                        xCoordsIsoBoundaryR = xB;
                        cellCBR = cellCBA;
                        
                    end
                    connecBound = cat(1,connecBound,connR);
                    xCoordsIsoBoundary = cat(3,xCoordsIsoBoundary,xCoordsIsoBoundaryR);
                    cellCB = cat(1,cellCB,cellCBR);
                    
                end
                
            end
            
            obj.connec                = connecT;
            obj.xCoordsIso            = xCoordsIso;
            obj.cellContainingSubcell = cellC;
            obj.xCoordsIsoBoundary = xCoordsIsoBoundary;
            obj.connecB = connecBound;
            obj.cellContainingSubCellBoundary = cellCB;
        end
        
        function [lNode,cNode] = obtainLongestAndClosestNode(obj,connec)
            xNodes  = obj.cutCoordComputer.coord;
            xNode1 = xNodes(connec(:,1),:);
            xNode2 = xNodes(connec(:,2),:);
            xNode3 = xNodes(connec(:,3),:);
            xNode4 = xNodes(connec(:,4),:);
            
            x12 = xNode2 - xNode1;
            x13 = xNode3 - xNode1;
            x14 = xNode4 - xNode1;
            
            nx12 = obj.computeNorm(x12);
            nx13 = obj.computeNorm(x13);
            nx14 = obj.computeNorm(x14);

            [~,indMax] = max([nx12,nx13,nx14],[],2);
            [~,indMin] = min([nx12,nx13,nx14],[],2);
            
            pos = [2;3;4];
            lNode = pos(indMax);
            cNode = pos(indMin);
        end
        
        function connec = obtainFirstTriangleConnectivity(obj,longestNode,closestNode)            
        %    connec = [ones(length(longestNode),1),longestNode,closestNode];
           connec = [ones(length(longestNode),1),2*ones(length(longestNode),1),4*ones(length(longestNode),1)];
        end
        
        function connec = obtainSecondTriangleConnectivity(obj,longestNode,closestNode,connec1) 
%             posT    = obj.computePosT(connec1);
%             [~,otherNode] = max(~posT,[],2);
%             connec = [ones(length(longestNode),1),otherNode,longestNode];
            connec = [ones(length(longestNode),1),3*ones(length(longestNode),1),4*ones(length(longestNode),1)];
        end
        
        function posT = computePosT(obj,connec1)
            posT = false(size(connec1,1),4);
            ind1 = sub2ind(size(posT),[1:size(connec1,1)]',connec1(:,1));
            ind2 = sub2ind(size(posT),[1:size(connec1,1)]',connec1(:,2));
            ind3 = sub2ind(size(posT),[1:size(connec1,1)]',connec1(:,3));
            posT(ind1) = true;
            posT(ind2) = true;
            posT(ind3) = true;
        end
        

      
   
        
        function n = computeNorm(obj,x)
            n = sqrt(x(:,1).^2 + x(:,2).^2 + x(:,3).^2);
        end
        
        function computeMesh(obj)
            sM.connec = obj.connec;
            sM.coord  = obj.coord;
            sM.kFace  = obj.backgroundMesh.kFace;
            obj.mesh = Mesh.create(sM);
        end
        
        function computeBoundaryMesh(obj)
            s.coord  = obj.cutCoordComputer.allValues;
            s.connec = obj.connecB;%obj.cutPointsInElemComputer.edgeCutPointInElem;
            s.kFace  = obj.backgroundMesh.kFace -1;
            obj.boundaryMesh = Mesh.create(s);
        end
        
%         function computeBoundaryXCoordsIso(obj)
%             obj.xCoordsIsoBoundary = obj.cutPointsInElemComputer.xCutInElem;
%         end

%         function computeBoundaryCellContainingSubCell(obj)
%             obj.cellContainingSubCellBoundary = obj.cutCells;
%         end
    end
    
end