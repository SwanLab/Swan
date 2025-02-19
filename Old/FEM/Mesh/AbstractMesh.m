classdef AbstractMesh < handle
    

    properties (GetAccess = public, SetAccess = protected)
        coord
        connec
        
        nelem
        ndim
        
        type
        
        coordElem
        interpolation
        
        edges
    end
    
    properties (Access = private)
       xFE
    end
    
    methods (Access = public)
       
       function xV = computeBaricenter(obj)
            xV = obj.xFE.computeValueInCenterElement();
       end
       
       function xGauss = computeXgauss(obj,xV)
            xGauss = obj.xFE.interpolateFunction(xV);
       end
       
       function dvolume = computeDvolume(obj,quad)
            s.mesh = obj;
            g = Geometry.create(s);
            g.computeGeometry(quad,obj.interpolation);
            dvolume = g.dvolu;
            dvolume = dvolume';
       end
       
       function q = computeElementQuality(obj)
            quad = Quadrature.set(obj.type);
            quad.computeQuadrature('CONSTANT');
            volume = obj.computeDvolume(quad); 
            L(1,:) = obj.computeSquarePerimeter();
            q = 4*sqrt(3)*volume./L;
       end
       
       function v = computeVolume(obj)
            quad = Quadrature.set(obj.type);
            quad.computeQuadrature('CONSTANT');
            v = obj.computeDvolume(quad);
            v = sum(v(:));
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
            obj.interpolation = Interpolation.create(obj.type,'LINEAR');
        end
        
        function computeElementCoordinates(obj)
            obj.computeCoordFEfunction();
            obj.coordElem = obj.xFE.fValues;
        end
        
    end
    
    methods (Access = private)
        
        function computeCoordFEfunction(obj)
            s.mesh    = obj;
            s.fValues = obj.coord;
            s.order   = 'P1';
            obj.xFE = LagrangianFunction(s);
        end
        
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