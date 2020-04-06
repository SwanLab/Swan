classdef AbstractMesh < handle
    
%     properties (Access = public)
%         unfittedType        
%     end

    
    
    properties (GetAccess = public, SetAccess = protected)
        coord
        connec
        
        nelem
        ndim
        
        geometryType
        
        coordElem
        interpolation
        
        edges
    end
    
    methods (Access = public)
       
       function xGauss = computeXgauss(obj,quad)
            xV     = quad.posgp;
            fElem  = obj.coordElem;
            fElem   = permute(fElem,[3 2 1]);                      
            xGauss = obj.interpolateFunction(xV,fElem);
       end      
       
       function fxV = interpolateFunction(obj,xV,func)
            intF   = obj.interpolation;           
            intF.computeShapeDeriv(xV);
            shapes = intF.shape;
            nNode  = size(shapes,1);
            nGaus  = size(shapes,2); 
            nElem  = size(func,1);                        
            nDime  = size(func,2);
            shapes = permute(shapes,[3 1 2]);            
            fxV = zeros(nElem,nDime,nGaus);
            for kNode = 1:nNode
                shapeKJ = shapes(:,kNode,:);
                fKJ(:,:,1) = func(:,:,kNode);
                f = bsxfun(@times,shapeKJ,fKJ);
                fxV = fxV + f;
            end
            fxV = permute(fxV,[2 3 1]);           
       end
       
       function dvolume = computeDvolume(obj,quad)
            s.mesh = obj;
            g = Geometry.create(s);
            g.computeGeometry(quad,obj.interpolation);
            dvolume = g.dvolu;
            dvolume = dvolume';
       end  
       
       function q = computeElementQuality(obj)
            quad = Quadrature.set(obj.geometryType);
            quad.computeQuadrature('CONSTANT');
            volume = obj.computeDvolume(quad);
            L(1,:) = obj.computeSquarePerimeter();
            q = 4*sqrt(3)*volume./L;           
       end
       
        function computeEdges(obj)
            s.nodesByElem = obj.connec;
            edge = EdgesConnectivitiesComputer(s);
            edge.compute();            
            obj.edges = edge;
        end       
                
    end
    
    methods (Access = protected)
        
        function createInterpolation(obj)
            obj.interpolation = Interpolation.create(obj,'LINEAR');            
        end        
        
        function computeElementCoordinates(obj)
            nNode  = obj.nnode;
            nDime  = obj.ndim;
            coordE = zeros(nNode,nDime,obj.nelem);
            coords = obj.coord';
            for inode = 1:nNode
                nodes = obj.connec(:,inode);
                coordNodes = coords(:,nodes);
                coordE(inode,:,:) = coordNodes;
            end
            obj.coordElem = coordE;
        end        
        
        
    end
    
    methods (Access = private)
        
        function L = computeSquarePerimeter(obj)
            obj.computeEdges();
            nElem = size(obj.connec,1);
            L = zeros(nElem,1);
            for iedge = 1:obj.edges.nEdgeByElem
                edge = obj.edges.edgesInElem(:,iedge);
                nodesEdge = obj.edges.nodesInEdges(edge,:);
                for idim = 1:obj.ndim
                    xA = obj.coord(nodesEdge(:,1),idim);
                    xB = obj.coord(nodesEdge(:,2),idim);
                    L = L + (xA - xB).^2;
                end
            end
        end
        
    end
    
end