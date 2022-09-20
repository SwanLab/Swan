classdef TestingCorrectors < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
       mesh
       orientation       
       singularityCoord
       boundaryPointCoord
       pathVertexes
       cellRight
       cellLeft
    end
    
    methods (Access = public)
        
        function obj = TestingCorrectors()
            obj.createMesh();
            obj.createOrientation();     
            obj.plotOrientationVector();
            obj.computeSingularities();
            obj.createBoundaryPoint();
            obj.computePathToBoundary();
            obj.createLeftRightPathElements();   
            obj.createSymmetricMapCond();
        end
        
    end
    
    methods (Access = private)
        
        function createMesh(obj)
            x = [66.90 128.89 115.25 76.73 26.84 157.24 ...
                168.58 141.74 98.65 45.74 3.78 2.65 31.37 ...
                83.53 137.21 174.63 147.03 99.03 38.55];
            y = [89.20 89.58 120.20 130.40 119.44 113.39 ...
                153.08 154.59 160.26 164.80 152.32 220.74 ...
                207.13 202.977 190.88 192.77 239.64 241.90 242.28];
            
         %   x = [0 1 1 0 0.5];
          %  y = [0 0 1 1 0.5];
            
            s.coord(:,1) = x;
            s.coord(:,2) = y;                        
            s.connec = delaunay(s.coord);
            m = Mesh(s);
            m.plot();
            obj.mesh = m;
        end
        
        function createOrientation(obj)
            alpha = pi/180*[190 140 300 150 185 ... 
                            120 75  45  150 130 ...
                            160+180 0   190  5  -5   ...
                            170 0   190 190];
            a(:,1) = cos(alpha);
            a(:,2) = sin(alpha);
            obj.orientation = a;
        end
        
        function plotOrientationVector(obj)
            %figure()
            a = obj.orientation;
            x = obj.mesh.coord(:,1);
            y = obj.mesh.coord(:,2);
            ax = a(:,1);
            ay = a(:,2);
            quiver(x,y,ax,ay);
        end        
        
        function computeSingularities(obj)
            s.mesh        = obj.mesh;
            s.orientation = obj.orientation;
            sF = SingularitiesFinder(s);
            isS = sF.computeSingularElements();
            sF.plot();
            coordB = obj.mesh.computeBaricenter();
            coordB = transpose(coordB);
            obj.singularityCoord = coordB(isS,:);
            obj.singularityCoord(:,2) = obj.singularityCoord(:,2) +0.1;
        end
        
        function createBoundaryPoint(obj)
            obj.boundaryPointCoord = [0 150];
        end
        
        function computePathToBoundary(obj)
            s.mesh = obj.mesh;
            s.singularityCoord   = obj.singularityCoord;
            s.boundaryPointCoord = obj.boundaryPointCoord;
            p = PathVertexesToBoundaryComputer(s);
            v = p.compute(); 
            obj.pathVertexes = v;
        end
        
        function createLeftRightPathElements(obj)
            s.pathVertexes = obj.pathVertexes;
            s.mesh         = obj.mesh;
            l = LeftRightCellsOfPathToBoundaryComputer(s);
            [cR,cL] = l.compute();   
            l.plot();            
            obj.cellLeft  = cL;
            obj.cellRight = cR;
        end         
        
        function createSymmetricMapCond(obj)   
            s.mesh        = obj.mesh.createDiscontinousMesh();
            s.orientation = obj.createDiscontinousField(obj.orientation);
            c = CoherentOrientationSelector(s);
            isC = c.isOrientationCoherent();
            
            isR = obj.cellRight;
            isL = obj.cellLeft;
            obj.mesh.connec(isR,:)
            
            v = obj.pathVertexes;
            
            
            for iNode = 1:3
                isC(isR,:)
                isC(isL,:)
            end
        end              
        
        function fD = createDiscontinousField(obj,fValues)
            s.connec = obj.mesh.connec;
            s.type   = obj.mesh.type;
            s.fNodes = fValues;
            f = FeFunction(s);            
            fD = f.computeDiscontinousField();
        end        
        
        
    end
    
end