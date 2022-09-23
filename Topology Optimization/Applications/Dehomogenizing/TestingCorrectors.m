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
       isCellLeft
       isCellRight
       isVertexInCell
       areVertexCoherent
       referenceCells
       symCond
       correctorValue
    end
    
    methods (Access = public)
        
        function obj = TestingCorrectors()
           obj.createMesh();            
           obj.createOrientation();                 
   %         obj.createBenchmarkMesh();            
   %         obj.createBenchmarkOrientation();            
            obj.plotOrientationVector();
            obj.computeSingularities();
            obj.createBoundaryPoint();
            obj.computePathToBoundary();
            obj.createLeftRightPathElements();   
            obj.computeCoherentOrientation();
            obj.createCorrector();
            obj.solveProjection();
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
        
       function createBenchmarkMesh(obj)
           h = 0.03;
           xmin = 0.50;
           xmax = 2.0;
           ymin = 0.25;
           ymax = 1.75;
           xv = xmin:h:xmax;
           yv = ymin:h:ymax;
           [X,Y] = meshgrid(xv,yv);
           s.coord(:,1) = X(:);
           s.coord(:,2) = Y(:);
           s.connec = delaunay(s.coord); 
%              [F,V] = mesh2tri(X,Y,zeros(size(X)),'x');
%              s.coord  = V(:,1:2);
%                s.connec = F;
            m = Mesh(s);            
            obj.mesh = m;
        end        
        
        function createBenchmarkOrientation(obj)
            s1 = 0.32;
            s2 = -0.8;
            x1 = obj.mesh.coord(:,1);
            x2 = obj.mesh.coord(:,2);         
            v(:,1) = cos(pi*(x1 + s1*x2));
            v(:,2) = cos(pi*(x2 + s2*x1));
            beta = atan2(v(:,2),v(:,1));
            alpha = beta/2;
            obj.orientation(:,1) = cos(alpha);
            obj.orientation(:,2) = sin(alpha);
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
            sCoord =  coordB(isS,:);
            obj.singularityCoord = sCoord(1,:);
            
            obj.singularityCoord(:,2) = sCoord(:,2)+0.1;
        end

        function createBoundaryPoint(obj)
            obj.boundaryPointCoord = [0 150];
       %     obj.boundaryPointCoord = [0.5 1.75];
        %   obj.boundaryPointCoord(:,1) = 1.75;            
        %   obj.boundaryPointCoord(:,2) = obj.singularityCoord(:,2);

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
            obj.isCellLeft  = cL;
            obj.isCellRight = cR;
        end                     
        
        function computeCoherentOrientation(obj)
            s.mesh        = obj.mesh.createDiscontinousMesh();
            s.orientation = obj.createDiscontinousField(obj.orientation);
            c = CoherentOrientationSelector(s);
            aC = c.isOrientationCoherent();
            obj.areVertexCoherent = aC;
        end        
        
       
        function createCorrector(obj)   
            s.isCellLeft        = obj.isCellLeft;
            s.isCellRight       = obj.isCellRight;
            s.areVertexCoherent = obj.areVertexCoherent;
            s.pathVertexes      = obj.pathVertexes;
            s.mesh              = obj.mesh;
            c = CorrectorComputer(s);
            cV = c.compute();
            c.plot()                        
            obj.correctorValue = cV;
        end
        
        function fGauss = createCorrectorGradient(obj)
            cV = obj.correctorValue;
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature('LINEAR');  
            
            int = Interpolation.create(obj.mesh,obj.mesh.interpolation.order);
            int.computeShapeDeriv(q.posgp); 
            dN = int.deriv;
            fGauss = zeros(obj.mesh.nelem,q.ngaus);
            for igaus = 1:q.ngaus
                for iNode = 1:size(dN,2)
                    grad = dN(igaus,iNode)*cV(:,iNode);
                    fGauss(:,igaus) = fGauss(:,igaus) + grad;
                end
            end            
        end
        
        function solveProjection(obj)

            s.mesh     = obj.mesh;
            s.fGauss   = obj.createCorrectorGradient();
            s.fValues  = obj.orientation';            
            
            m = MinimumDiscGradFieldWithVectorInL2(s);
            f = m.solve();
        end
        
        function createSymmetricMapCond(obj)   
            s.meshCont  = obj.mesh;
            s.meshDisc  = obj.mesh.createDiscontinousMesh();
            s.fieldDisc =  obj.createDiscontinousField(obj.orientation);
            s = SymmetricContMapCondition(s);            
            sC = s.computeCondition();
            obj.symCond = sC;            
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