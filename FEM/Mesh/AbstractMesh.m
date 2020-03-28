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
    end
    
    methods (Access = public)
       
       function xGauss = computeXgauss(obj,quad)
            xV = quad.posgp;
            xGauss = obj.interpolateFunction(xV,obj.coordElem);
       end      
       
       function fxV = interpolateFunction(obj,xV,func)
            obj.interpolation.computeShapeDeriv(xV);           
            nNode  = obj.nnode;
            nDime  = size(func,2);
            shapes = obj.interpolation.shape;
            nGaus  = size(xV,2);
            fxV = zeros(nGaus,nDime,obj.nelem);
            for kNode = 1:nNode
                shapeKJ(:,1) = shapes(kNode,:)';
                fKJ = func(kNode,:,:);
                f = bsxfun(@times,shapeKJ,fKJ);
                fxV = fxV + f;
            end
            fxV = permute(fxV,[2 1 3]);           
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
            nodes = obj.connec;
            s.nodesByElem = nodes;
            edges = EdgesConnectivitiesComputer(s);
            edges.compute();
            nElem = size(nodes,1);
            L = zeros(nElem,1);
            for iedge = 1:edges.nEdgeByElem
                edge = edges.edgesInElem(:,iedge);
                nodesEdge = edges.nodesInEdges(edge,:);
                for idim = 1:obj.ndim
                    xA = obj.coord(nodesEdge(:,1),idim);
                    xB = obj.coord(nodesEdge(:,2),idim);
                    L = L + (xA - xB).^2;
                end
            end
        end
       
    end
    
end