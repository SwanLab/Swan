classdef Mesh < handle
    
    properties (GetAccess = public, SetAccess = private)
        nnode
        npnod
        type
        kFace    
    
        coord
        connec
        
        nelem
        ndim
        
        
        coordElem
        interpolation
        
        edges
    end
    
    properties (Access = private)
       xFE 
    end
    
    methods (Access = public)
        
        function obj = Mesh(cParams)
            obj.init(cParams);
            obj.computeDescriptorParams();            
            obj.computeType();
            obj.createInterpolation();
            obj.computeElementCoordinates();
        end
        
        function obj = createFromFile(obj,cParams)
            testName = cParams.testName;
            [coordV, connecV] = obj.readCoordConnec(testName);
            s.coord  = coordV(:,2:end-1);
            s.connec = connecV(:,2:end);
            obj = Mesh(s);
        end
  
        function plot(obj)
            s.mesh = obj;
            s = SettingsMeshPlotter(s);
            mP = MeshPlotter(s);
            mP.plot();
        end        
        
        function hMin = computeMinCellSize(obj)
            x1(:,1) = obj.coord(obj.connec(:,1),1);
            x1(:,2) = obj.coord(obj.connec(:,1),2);
            x2(:,1) = obj.coord(obj.connec(:,2),1);
            x2(:,2) = obj.coord(obj.connec(:,2),2);
            x3(:,1) = obj.coord(obj.connec(:,3),1);
            x3(:,2) = obj.coord(obj.connec(:,3),2);            
            x1x2 = (x2-x1);
            x2x3 = (x3-x2);
            x1x3 = (x1-x3);
            n12 = sqrt(x1x2(:,1).^2 + x1x2(:,2).^2);
            n23 = sqrt(x2x3(:,1).^2 + x2x3(:,2).^2);
            n13 = sqrt(x1x3(:,1).^2 + x1x3(:,2).^2);
            hs = min([n12,n23,n13],[],2);
            hMin = min(hs);            
        end
        
        function hMean = computeMeanCellSize(obj)
            x1(:,1) = obj.coord(obj.connec(:,1),1);
            x1(:,2) = obj.coord(obj.connec(:,1),2);
            x2(:,1) = obj.coord(obj.connec(:,2),1);
            x2(:,2) = obj.coord(obj.connec(:,2),2);
            x3(:,1) = obj.coord(obj.connec(:,3),1);
            x3(:,2) = obj.coord(obj.connec(:,3),2);            
            x1x2 = (x2-x1);
            x2x3 = (x3-x2);
            x1x3 = (x1-x3);
            n12 = sqrt(x1x2(:,1).^2 + x1x2(:,2).^2);
            n23 = sqrt(x2x3(:,1).^2 + x2x3(:,2).^2);
            n13 = sqrt(x1x3(:,1).^2 + x1x3(:,2).^2);
            hs = max([n12,n23,n13],[],2);
            hMean = max(hs);
        end
        
        function L = computeCharacteristicLength(obj)
            xmin = min(obj.coord);
            xmax = max(obj.coord);
            L = norm(xmax-xmin);%/2;
        end
        
        function changeCoordinates(obj,newCoords)
            obj.coord = newCoords;
        end
        
        function setCoord(obj,newCoord)
            obj.coord = newCoord;
        end
        
        function setConnec(obj,newConnec)
            obj.connec = newConnec;
        end
        
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
       
        function computeDescriptorParams(obj)
            obj.npnod = size(obj.coord,1);
            obj.ndim  = size(obj.coord,2);
            obj.nelem = size(obj.connec,1);
            obj.nnode = size(obj.connec,2);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            s = SettingsMesh(cParams);            
            obj.coord  = s.coord;
            obj.connec = s.connec;
            obj.type   = s.type;
            obj.kFace  = s.kFace;
        end        
        
        function computeType(obj)
            s.ndim = obj.ndim;
            s.nnode = obj.nnode;
            s.kFace = obj.kFace;
            t = MeshTypeComputer(s);
            obj.type = t.compute();
        end
        
        function createInterpolation(obj)
            obj.interpolation = Interpolation.create(obj,'LINEAR');            
        end        
        
        function computeElementCoordinates(obj)
            obj.computeCoordFEfunction();
            obj.coordElem = obj.xFE.fElem;
        end             
        
        function computeCoordFEfunction(obj)
            s.mesh   = obj;
            s.fNodes = obj.coord;
            obj.xFE = FeFunction(s);             
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
    
    methods (Access = private, Static)
       
        function [coord, connec] = readCoordConnec(testName)
            run(testName)            
        end     
        
    end
    
end

