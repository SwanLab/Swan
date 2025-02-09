classdef NextVertexToBoundaryComputer < handle
       
    properties (Access = private)
        currentVertex
        trialVertexes
        lineVector
        currentVertexCoord
        trialVertexesCoord
        boundaryVertexCoord         
    end
    
    methods (Access = public)
        
        function obj = NextVertexToBoundaryComputer(cParams)
            obj.init(cParams)            
        end
        
        function iD = compute(obj)
            merit   = obj.computeMeritFunction();
            [~,iD]  = min(merit); 
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.currentVertexCoord  = cParams.currentVertexCoord;
            obj.trialVertexesCoord  = cParams.trialVertexesCoord;
            obj.boundaryVertexCoord = cParams.boundaryVertexCoord;
            obj.lineVector          = cParams.lineVector;            
        end
        
        function merit = computeMeritFunction(obj)
            distToLine  = obj.computeDistanceToLine();
            distToBound = obj.computeDistanceToBoundary();
            indicator   = obj.computeIsNewPathAdvancing();
            alpha = 0.5;
            merit = alpha*distToLine + (1-alpha)*distToBound + indicator;                                          
        end
        
        function vOrtNorm = computeDistanceToLine(obj)
            v    = obj.computeVertexToBoundaryDirection();
            vOrt = obj.computeOrthogonalLineProjection(v);
            vOrtNorm = obj.computeVectorNorm(vOrt);
        end
        
        function vNorm = computeDistanceToBoundary(obj)
            v   = obj.computeVertexToBoundaryDirection();
            vNorm = obj.computeVectorNorm(v);
        end        
        
        function dir = computeVertexToBoundaryDirection(obj)            
            coordTrial = obj.trialVertexesCoord;
            coordBound = obj.boundaryVertexCoord;
            dir(:,1) = coordBound(:,1)  - coordTrial(:,1);
            dir(:,2) = coordBound(:,2)  - coordTrial(:,2);
        end              
        
        function indicator = computeIsNewPathAdvancing(obj)
            edge        = obj.computeNewEdgeDirection();
            isAdvancing = obj.isEdgeOrientedWithPathVector(edge);
            indicator   = zeros(size(isAdvancing));
            indicator(~isAdvancing) = inf;
        end

        function uOrt = computeOrthogonalLineProjection(obj,u)
            vectorLine = obj.lineVector;
            uProj = obj.computeProjection(u,vectorLine);
            uOrt  = u - uProj;
        end
        
        function uProj = computeProjection(obj,u,vectorLine)
            w     = obj.computeUnitVector(vectorLine);
            uProj = obj.scalarProduct(u,w)*w;
        end   
        
        function dir = computeNewEdgeDirection(obj)
            OA = obj.currentVertexCoord;
            OB = obj.trialVertexesCoord;
            AB(:,1) = OB(:,1)  - OA(1);
            AB(:,2) = OB(:,2)  - OA(2);
            dir = AB;
        end        
        
        function itIs = isEdgeOrientedWithPathVector(obj,edgeDir)
            a = obj.lineVector;
            b = edgeDir;
            itIs = obj.areVectorsEquallyOriented(a,b);
        end
        
        function theyAre = areVectorsEquallyOriented(obj,a,b)
            ab  = obj.scalarProduct(b,a);
            theyAre = ab > 0;
        end
        
        function vNorm = computeVectorNorm(obj,v)
            vv = obj.scalarProduct(v,v);
            vNorm  = sqrt(vv);
        end        
        
        function u = computeUnitVector(obj,AB)
            nAB = obj.computeVectorNorm(AB);
            u(:,1) = AB(:,1)./nAB;
            u(:,2) = AB(:,2)./nAB;
        end
                
    end
    
    methods (Access = private, Static)
        
        function ab = scalarProduct(a,b)
            ab = a(:,1).*b(:,1) + a(:,2).*b(:,2);
        end
                
    end    
    
end