% classdef IntersectionCoordComputer < TotalCoordinatesCalculator
% 
%     properties (Access = public)
%         c
%         theta
%         nodes
%         vertCoord
%         boundCoord
%         div
%         totalCoord
%     end
% 
%     properties (Access = private)
%         vA
%         vB
%     end
% 
%     methods (Access = public)
% 
%         function obj = IntersectionCoordComputer(cParams)
%             obj.init(cParams);
%             obj.initBoundary();
%             obj.obtainPrincipalVectors();
%             obj.computeIntersections();
%         end
% 
%     end
% 
%     methods (Access = private)
% 
%         function obtainPrincipalVectors(obj)
%             A = obj.vertCoord(1,:);
%             B = obj.vertCoord(2,:);
%             C = obj.vertCoord(3,:);
%             obj.vA = IntersectionCoordComputer.computeUnitaryVector(A,B);
%             obj.vB = IntersectionCoordComputer.computeUnitaryVector(B,C);
%         end
% 
%         function computeIntersections(obj)
%             nodesX = obj.div(1)-1;
%             nodesY = obj.div(2)-1;
%             intNode = obj.nodes.bound +1;
%             for jNodes = 1:nodesY
%                 pB = obj.boundCoord(obj.nodes.bound+1-jNodes,:);
%                 for iNodes = 1:nodesX
%                     pA = obj.boundCoord(obj.nodes.vert+iNodes,:);
%                     if obj.vA(1) == 0
%                         [x,y] = IntersectionCoordComputer.computeVerticalIntersection(pA,pB,obj.vB);
%                     elseif obj.vB(1) == 0
%                         [x,y] = IntersectionCoordComputer.computeVerticalIntersection(pB,pA,obj.vA);
%                     else
%                         [x,y] = IntersectionCoordComputer.computeGeneralIntersection(pA,pB,obj.vA,obj.vB);
%                     end
%                     obj.totalCoord(intNode,:) = [x y];
%                     intNode = intNode+1;
%                 end
%             end
%         end
% 
%     end
% 
%     methods (Static)
% 
%         function v_norm = computeUnitaryVector(A,B)
%             v = B-A;
%             m = norm(v);
%             v_norm = v/m;
%         end
% 
%         function [x,y] = computeVerticalIntersection(p1,p2,v)
%             x = p2(1);
%             y = (x-p1(1))*v(2)/v(1)+p1(2);
%         end
% 
%         function [x,y] = computeGeneralIntersection(p1,p2,v1,v2)
%             x = (p2(2)-p1(2)+p1(1)*v2(2)/v2(1)-p2(1)*v1(2)/v1(1))/(v2(2)/v2(1)-v1(2)/v1(1));
%             y = (x-p2(1))*v1(2)/v1(1)+p2(2);
%         end
% 
%     end
% 
% end

classdef IntersectionCoordComputer < TotalCoordinatesCalculator
    
    properties (Access = public)
        latticeVectors
        nodes
        vertCoord
        boundCoord
        div
        totalCoord
    end
    
    properties (Access = private)
        vA
        vB
    end
    
    methods (Access = public)
        
        function obj = IntersectionCoordComputer(cParams)
            obj.init(cParams);
            obj.initBoundary();
            obj.obtainPrincipalVectors();
            obj.computeIntersections();
        end
        
    end
    
    methods (Access = private)
        
        function obtainPrincipalVectors(obj)
            
            if ~isempty(obj.latticeVectors) && size(obj.latticeVectors, 1) == 2
                
                v1 = obj.latticeVectors(1,:);
                v2 = obj.latticeVectors(2,:);
                obj.vA = v1 / norm(v1);
                obj.vB = v2 / norm(v2);
            else
                
                A = obj.vertCoord(1,:);
                B = obj.vertCoord(2,:);
                C = obj.vertCoord(3,:);
                obj.vA = IntersectionCoordComputer.computeUnitaryVector(A,B);
                obj.vB = IntersectionCoordComputer.computeUnitaryVector(B,C);
            end
        end
        
        function computeIntersections(obj)
            nDiv1 = obj.div(1);
            nDiv2 = obj.div(2);
            nVert = obj.nodes.vert;  
            
            
            startL1 = nVert + 1;
            
            
            startL2 = startL1 + (nDiv1 - 1);
            
          
            
            intNode = obj.nodes.bound + 1;
            
           
            for jNodes = 1:nDiv2-1
               
                pB = obj.boundCoord(startL2 + jNodes - 1, :);
                
                for iNodes = 1:nDiv1-1
                    
                    pA = obj.boundCoord(startL1 + iNodes - 1, :);
                    
                    if abs(obj.vA(1)) < 1e-12
                        [x,y] = IntersectionCoordComputer.computeVerticalIntersection(pA, pB, obj.vB);
                    elseif abs(obj.vB(1)) < 1e-12
                        [x,y] = IntersectionCoordComputer.computeVerticalIntersection(pB, pA, obj.vA);
                    else
                        [x,y] = IntersectionCoordComputer.computeGeneralIntersection(pA, pB, obj.vA, obj.vB);
                    end
                    
                    obj.totalCoord(intNode,:) = [x, y];
                    intNode = intNode + 1;
                end
            end
        end
        
    end
    
    methods (Static)
        
        function v_norm = computeUnitaryVector(A,B)
            v = B - A;
            m = norm(v);
            if m > 1e-12
                v_norm = v / m;
            else
                v_norm = [0, 0];
            end
        end
        
        function [x,y] = computeVerticalIntersection(p1, p2, v)
            x = p2(1);
            if abs(v(1)) > 1e-12
                y = (x - p1(1)) * v(2)/v(1) + p1(2);
            else
                y = p1(2);
            end
        end
        
        function [x,y] = computeGeneralIntersection(p1, p2, v1, v2)
            denom = v2(2)/v2(1) - v1(2)/v1(1);
            if abs(denom) > 1e-12
                x = (p2(2) - p1(2) + p1(1)*v2(2)/v2(1) - p2(1)*v1(2)/v1(1)) / denom;
                y = (x - p2(1)) * v1(2)/v1(1) + p2(2);
            else
                
                x = (p1(1) + p2(1)) / 2;
                y = (p1(2) + p2(2)) / 2;
            end
        end
        
    end
    
end