classdef QuadrilaterSubTriangulator < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        subTriangulation
    end
    
    properties (Access = private)
        localIncidenceMatrix
        coord
        connec
        xLocalVertex
    end
    
    methods (Access = public)
        
        function obj = QuadrilaterSubTriangulator(cParams)
            obj.init(cParams)
        end
        
        function [connR,xCoordsIsoBoundaryR] = compute(obj)
            
            for is = 1:2
                subTriang = obj.subTriangulation(:,:,is);
                connecL             = obj.computeSubtriangulation(subTriang);
                conne               = obj.divideInSubTriangle(obj.connec,connecL);

                xB = permute(obj.xLocalVertex,[3 2 1]);
                xCoordsIsoBoundaryS = obj.divideInSubTriangle(xB,connecL);
                
                xCoordsIsoBoundaryS = permute(xCoordsIsoBoundaryS,[4 1 2 3]);
                
                conn1(:,:,is) = conne(:,:,1);
                conn2(:,:,is) = conne(:,:,2);
                %connRi(:,:,is) = cat(1,conn1,conn2);
                
                
                nodesSubCases(:,:,is,:) = permute(conne,[3 2 1]);
                
                
                xCoordsIsoBoundary1(:,:,:,is) = xCoordsIsoBoundaryS(:,:,:,1);
                xCoordsIsoBoundary2(:,:,:,is) = xCoordsIsoBoundaryS(:,:,:,2);
                
                
                %  nodesSubCases(:,:,1,:) =
            end
            

            
            
            s.nodesSubCases = nodesSubCases;
            s.coord = obj.coord;
            b = BestSubCellCaseSelector(s);
            imax = b.compute();
            
            conR1 = zeros(size(conn1(:,:,1)));
            conR2 = zeros(size(conn1(:,:,1)));
            xCoordsIsoBoundaryR1 = zeros(size(xCoordsIsoBoundary1(:,:,:,1)));
            xCoordsIsoBoundaryR2 = zeros(size(xCoordsIsoBoundary1(:,:,:,1)));
            for is = 1:2
                isIndex = imax == is;
                con1 = conn1(:,:,is);
                con2 = conn2(:,:,is);
                
                conR1(isIndex,:) = con1(isIndex,:);
                conR2(isIndex,:) = con2(isIndex,:);
                
                xR1 = xCoordsIsoBoundary1(:,:,:,is);
                xR2 = xCoordsIsoBoundary2(:,:,:,is);
                xCoordsIsoBoundaryR1(:,isIndex,:) = xR1(:,isIndex,:);
                xCoordsIsoBoundaryR2(:,isIndex,:) = xR2(:,isIndex,:);
                
            end
            %   connR = cat(1,conR1,conR2);
            connR = cat(1,conne(:,:,1),conne(:,:,2));
            
            
            xCoordsIsoBoundaryRN = cat(2,xCoordsIsoBoundaryR1,xCoordsIsoBoundaryR2);
            xCoordsIsoBoundaryR = permute(xCoordsIsoBoundaryRN,[1 3 2]);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.localIncidenceMatrix = cParams.localIncidenceMatrix;
            obj.coord                = cParams.coord;
            obj.connec               = cParams.connec;
            obj.xLocalVertex         = cParams.xLocalVertex;
            obj.subTriangulation(:,:,1) = [1 2 3; 3 4 1];
            obj.subTriangulation(:,:,2) = [2 3 4; 4 1 2];
        end
        
        function [xCoordsIsoBoundaryR] = divideInSubTriangle(obj,xB,connecL)
            
            
     %       connB = obj.localToGlobal(connecL,connec);
            
            nElem = size(connecL,1);
            xCoordsIsoBoundaryR = zeros(nElem,3,2,size(xB,3));
            for idim = 1:size(xB,3)
                xBi = (xB(:,:,idim));
                xC = obj.localToGlobal(connecL,xBi);
                xCoordsIsoBoundaryR(:,:,:,idim) = xC;
            end
        end
        
        
        
        function [connecG] = localToGlobal(obj,connecL,dataToTriang)
            
            for itriang = 1:2
                for i = 1:3
                    nodei1 = connecL(:,i,itriang);
                    index1 = sub2ind(size(dataToTriang),[1:size(connecL,1)]',nodei1);
                    connecG(:,i,itriang) = dataToTriang(index1);
                end
            end
            
        end
        
        function nodes = obtainOrderedNodes(obj)
            incMatrix = obj.localIncidenceMatrix;
            nElem = size(incMatrix,1);
            incMatrix = permute(incMatrix,[1 3 2]);
            edge1 = incMatrix(:,:,1);
            edge2 = incMatrix(:,:,2);
            edge3 = incMatrix(:,:,3);
            edge4 = incMatrix(:,:,4);
            areIntersected(:,1) = obj.areEdgesIntersected(edge1,edge2);
            areIntersected(:,2) = obj.areEdgesIntersected(edge1,edge3);
            areIntersected(:,3) = obj.areEdgesIntersected(edge1,edge4);
            
            integerCases = obj.computeIntegerCases(areIntersected);
            
            nodes = zeros(nElem,4);
            cases = [3 5 6];
            nodesCases = [1 2 4 3; 1 2 3 4; 1 3 2 4];
            for icase = 1:3
                isCase = integerCases == cases(icase);
                for inode = 1:4
                    nodes(isCase,inode) = nodesCases(icase,inode);
                end
            end
            
        end
        
        function connec = computeSubtriangulation(obj,subTriang)
            nodes = obj.obtainOrderedNodes();
            connec = zeros(size(nodes,1),3,2);
            for inode = 1:3
                for itriangle = 1:2
                    connec(:,inode,itriangle) = nodes(:,subTriang(itriangle,inode));
                end
            end
        end
        
    end
    
    methods (Access = private, Static)
        

        
        function [areInt] = areEdgesIntersected(edge1,edge2)
            isNode11 = edge1(:,1) == edge2(:,1);
            isNode12 = edge1(:,1) == edge2(:,2);
            isNode21 = edge1(:,2) == edge2(:,1);
            isNode22 = edge1(:,2) == edge2(:,2);
            areInt = any([isNode11 isNode12 isNode21 isNode22],2);
        end
        
        function integerCases = computeIntegerCases(isNeigbour)
            nodes = isNeigbour;
            nnode = size(nodes,2);
            nodePos = (1:nnode) - 1;
            pow2vector = 2.^(nodePos);
            d = pow2vector*nodes';
            integerCases = d;
        end
        
        
        
    end
    
end