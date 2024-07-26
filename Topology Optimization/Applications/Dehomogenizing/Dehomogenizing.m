classdef Dehomogenizing < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        backgroundMesh
        boundaryMesh
        uMesh
        coord 
        theta
    end
    
    properties (Access = private)
        nx1
        nx2
        nCells
    end
    
    methods (Access = public)
        
        function obj = Dehomogenizing()
            obj.init();
            obj.createCoord();
            obj.createBackgroundMesh();
            obj.createOrientation();
            obj.createBoundaryMesh();
            obj.createUnfittedMesh();
            obj.plotOrientation();
            obj.plotStructure();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.nx1    = 25*2;
            obj.nx2    = 25;
            obj.nCells = 20;
        end
        
        function createCoord(obj)
            x1 = linspace(-1,1,obj.nx1);
            x2 = linspace(0,1,obj.nx2);
            x1T = repmat(x1,obj.nx2,1);
            x2T = repmat(x2',1,obj.nx1);
            obj.coord = [x1T(:),x2T(:)];
        end
        
        function createOrientation(obj)
            x2 = obj.coord(:,2);
            x1 = obj.coord(:,1);
            obj.theta = atan2(x2,x1);
        end
        
        function createBackgroundMesh(obj)
%             xy = obj.coord;
%             x1 = xy(:,1);
%             x2 = xy(:,2);
%             s.connec   = delaunay(x1,x2);
%             s.coord    = obj.coord;
%             obj.backgroundMesh = Mesh(s);
%             
            
            x1 = linspace(-1,1,obj.nx1);
            x2 = linspace(0,1,obj.nx2);
            x1min = min(x1);
            x1max = max(x1);
            x2min = min(x2);
            x2max = max(x2);
            [coordinates, nodes,nel,nnode] = MeshRectangular(x1max-x1min,x2max-x2min,obj.nx1,obj.nx2);
            s.coord(:,1) = coordinates(:,1)+x1min;
            s.coord(:,2) = coordinates(:,2)+x2min;
            s.connec = nodes;
            obj.backgroundMesh = Mesh.create(s);  
            obj.backgroundMesh.plot()
            obj.coord = s.coord;
            
            
%             x1 = linspace(-1,1,obj.nx1);
%             x2 = linspace(0,1,obj.nx2);
%             [xv,yv] = meshgrid(x1,x2); 
%              [F,V] = mesh2tri(xv,yv,zeros(size(xv)),'x');
%              s.coord  = V(:,1:2);
%              s.connec = F;
%              obj.backgroundMesh = Mesh(s);
%              obj.backgroundMesh.plot()
%              obj.coord = s.coord;
            
        end
        
        function plotOrientation(obj)
             figure()
             x = obj.backgroundMesh.coord(:,1);
             y = obj.backgroundMesh.coord(:,2);
             t  = obj.theta;
             ct = cos(t(:,1));
             st = sin(t(:,1));
             quiver(x,y,ct,st)
        end
        
        function plotStructure(obj)
            figure()
            obj.uMesh.plotStructureInColor('black');
        end
        
        function createBoundaryMesh(obj)
            sB.backgroundMesh = obj.backgroundMesh;
            sB.dimension = 1:3;
            sB.type = 'FromReactangularBox';
            bMc = BoundaryMeshCreator.create(sB);
            obj.boundaryMesh  = bMc.create();
        end
        
        function createUnfittedMesh(obj)
            ls = obj.createLevelSet();
            s.boundaryMesh   = obj.boundaryMesh;
            s.backgroundMesh = obj.backgroundMesh;
            obj.uMesh = UnfittedMesh(s);
            obj.uMesh.compute(ls)
        end
        
        function ls = createLevelSet(obj)
            s.coord = obj.coord;
            s.type   = 'periodicAndOriented';
            s.mesh   = obj.backgroundMesh;
            s.ndim   = 2;
            s.angle  = obj.theta;
            s.nCells = obj.nCells;
            s.cellLevelSetParams = obj.createLevelSetCellParams();
            levelSet = LevelSetCreator.create(s);
            ls = levelSet.getValue();
        end
        
        function s = createLevelSetCellParams(obj)
           sE = obj.createSuperEllipseParams();
           s.type   = 'smoothRectangle';
           s.widthH = sE.m1;
           s.widthV = sE.m2;
           s.pnorm  = sE.q;
           s.ndim   = 2;
        end
        
        function sE = createSuperEllipseParams(obj)
           s.coord = obj.coord;
           s.mMin  = 0.58;
           s.mMax  = 0.28;
           s.qMin  = 36;
           s.qMax  = 36;
           sE = SuperEllipseDistributionExample(s);
           sE.computeParameters();
        end
    
    end
    
end