classdef PathVertexesToBoundaryComputer < handle
    
    properties (Access = private)
       singularityCoord
       boundaryPointCoord               
       pathVector
       initialVertex
       boundaryNodes
       connectedVertex     
       pathVertexes       
       pathCells
    end

    properties (Access = private)
       mesh 
       singularElement
       isCoherent       
    end    
    
    methods (Access = public)
        
        function obj = PathVertexesToBoundaryComputer(cParams)
            obj.init(cParams)
        end
        
        function v = compute(obj)
            obj.computeInitialVertex();            
            obj.computeBoundaryPointCoord();            
            obj.computeSingularityCoord();            
            obj.createStraightPathVector();
            obj.computeBoundaryNodes();
            obj.computeMeshEdges();
            obj.computePath();
            v = obj.pathVertexes;          
        end

        function plot(obj)
            figure()
            hold on
            obj.plotMesh();
            obj.plotStraightPath();
            obj.plotInitialVertex();
            obj.plotVerticesPath();
        end        
        
    end
    
    methods (Access = private)
              
        
        function init(obj,cParams)
            obj.mesh               = cParams.mesh;
            obj.singularElement   = cParams.singularElement;
            obj.isCoherent        = cParams.isCoherent;
        end
        
        function sC = computeSingularElementBaricenter(obj)
             isS = obj.singularElement;
             sC  = obj.mesh.computeBaricenter();
             sC  = transpose(sC);
             sC  = sC(isS,:);              
        end

        function computeSingularityCoord(obj)
            obj.singularityCoord = obj.mesh.coord(obj.initialVertex,:);
        end

        function computeDistance(a,b)

        end

        function computeBoundaryPointCoord(obj)
            sC = obj.computeSingularElementBaricenter();
            % bM = obj.mesh.createBoundaryMesh();
            % for iB = 1:numel(bM)
            %    bmI = bM{iB};
            %    bC  = bmI.mesh.coord;
            %   dis = (bC(:,1) - sC(1)).^2 + (bC(:,2) - sC(2)).^2;
            %   % dis = (bC(:,1) - sC(1)).^2; %+ (bC(:,2) - sC(2)).^2;               
            %    [dM,im] = min(dis);
            %    nodes(iB) = bmI.globalConnec(im);
            %    dist(iB) = dM;
            % end

            nodes = boundary(obj.mesh.coord,1);
            bC = obj.mesh.coord(nodes,:);
            dist = (bC(:,1) - sC(1)).^2 + (bC(:,2) - sC(2)).^2;


            [~,iM] = min(dist);
            node = nodes(iM);
            cP = obj.mesh.coord(node,:);
       %   obj.boundaryPointCoord = [1.98    1.3300];
          % obj.boundaryPointCoord = [0.5    1.3900];
    %       obj.boundaryPointCoord = [2    0.55];
           obj.boundaryPointCoord = cP;
          
            %coordI = obj.mesh.coord(obj.initialVertex,:);
            %xmin = min(obj.mesh.coord(:,1));
            %obj.boundaryPointCoord = [xmin coordI(2)];
            
        end
                  

        function createStraightPathVector(obj)
            OA = obj.singularityCoord;            
            OB = obj.boundaryPointCoord;
            AB = OB - OA;
            obj.pathVector = AB ;
        end   

        function itIs = isCellAroundCoherent(obj,vertex)
            cells = obj.mesh.computeAllCellsOfVertex(vertex);
            isCV  = obj.isCoherent.getFvaluesDisc();
            isCo = squeeze(isCV(1,:,cells));
            itIs = all(isCo);
        end

        function computeInitialVertex(obj)
            isS = obj.singularElement;
            %  sN  = obj.mesh.connec(isS,:);
            %  xy  = obj.mesh.coord(sN,:);
            %  xyS = obj.boundaryPointCoord;
            %  xD  = xy(:,1) - xyS(:,1);
            %  yD  = xy(:,2) - xyS(:,2);
            %  dist = sqrt(xD.^2 + yD.^2);
            %  [~,iD] = min(dist);
            % 
            %  for idof = 1:obj.mesh.nnodeElem
            %     isC(idof) = obj.isCoherent.isDofContinous(isS,idof);
            %  end
            % iD = find(isC);

            for idof = 1:obj.mesh.nnodeElem
                iVertex = obj.mesh.connec(isS,idof);
                isCellCoherent = obj.isCellAroundCoherent(iVertex);
                notCoh(idof) = sum(~isCellCoherent);
            end
            
            isCV = obj.isCoherent.getFvaluesDisc();
            if all(isCV(1,:,isS))
                iD = find(notCoh == 1);
               % [~,iD] = min(notCoh);
            else
               [~,iD] = min(notCoh);
            end
            firstiD = iD;
          % iD = 3;
            cV     = obj.mesh.connec(isS,firstiD);
            obj.initialVertex = cV;
        end

        function computeSecondVertex(obj)
            
        end

        function computeBoundaryNodes(obj)
            coord  = obj.mesh.coord;
            nodesB = boundary(coord,1);            
            obj.boundaryNodes = nodesB;
        end
        
        function computeMeshEdges(obj)
            obj.mesh.computeEdges();            
        end

        function itIs = isInBoundary(obj,node)
            nodesB = obj.boundaryNodes;
            itIs = any(node == nodesB);
        end
        
        function itIs = itIsInFracture(obj,vertex)
           isAroundCellCoherent = obj.isCellAroundCoherent(vertex);
           itIs = any(~isAroundCellCoherent);
        end

        function d = distanceToBoundary(obj,vertex)
            bc = obj.boundaryPointCoord;
            cV = obj.mesh.coord(vertex,:);
            d = (bc(1) - cV(1))^2 + (bc(2) - cV(2))^2;
        end

        function d = distanceToSingularity(obj,vertex)
            sc = obj.computeSingularElementBaricenter();            
            cV = obj.mesh.coord(vertex,:);
            d = (sc(1) - cV(1))^2 + (sc(2) - cV(2))^2;
        end        

        function computePath(obj)            
            i = 1;
            vertex       = obj.initialVertex;
            pVertexes    = vertex;
            while ~obj.isInBoundary(vertex)
                otherVertex = obj.mesh.computeConnectedVertex(vertex);
                iD          = obj.computeOptimalVertex(vertex,otherVertex);  
                newVertex   = otherVertex(iD);  
                if i == 1
                   if obj.itIsInFracture(newVertex)
                    oVertex =  obj.mesh.computeConnectedVertex(obj.initialVertex);
                    itIsNewVertex = newVertex == oVertex;
                    candidates = oVertex(itIsNewVertex);
                    sC = obj.computeSingularElementBaricenter();                    
                    for iVer = 1:length(candidates)
                        cand = candidates(iVer);
                        itIsF(iVer) = obj.itIsInFracture(cand);
                        dB = obj.distanceToBoundary(cand);
                        dS = obj.distanceToSingularity(cand);
                   end 
                   [~,isClosest] = sort(dB);
                   end
                end
                i = i + 1;  
                vertex       = newVertex;                
                pVertexes(i) = vertex;              
            end
            obj.pathVertexes = pVertexes;
        end        
        
        function iD = computeOptimalVertex(obj,currentVertex,trialVertexes)             
            s.currentVertexCoord  = obj.mesh.coord(currentVertex,:);
            s.trialVertexesCoord  = obj.mesh.coord(trialVertexes,:);
            s.boundaryVertexCoord = obj.boundaryPointCoord;
            s.lineVector          = obj.pathVector;
            n = NextVertexToBoundaryComputer(s);
            iD = n.compute();
        end            
        
        function plotMesh(obj)
            obj.mesh.plot();
        end
        
        function plotStraightPath(obj)
            xP = [obj.singularityCoord(:,1),obj.boundaryPointCoord(:,1)];
            yP = [obj.singularityCoord(:,2),obj.boundaryPointCoord(:,2)];
            zP = ones(size(xP));
            plot3(xP,yP,zP,'+-')
        end
        
        function plotInitialVertex(obj)
            cV = obj.initialVertex;
            x = obj.mesh.coord(cV,1);
            y = obj.mesh.coord(cV,2);
            z = ones(size(x));
            plot3(x,y,z,'r+')
        end
        
        function plotVerticesPath(obj)
            cV = obj.pathVertexes;
            x = obj.mesh.coord(cV,1);
            y = obj.mesh.coord(cV,2);
            z = ones(size(x));
            plot3(x,y,z,'g-','LineWidth',5)
        end
        
    end
    
    
    
end