classdef Dehomogenizing < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        backgroundMesh
        boundaryMesh
        uMesh
        gMesh
        dMesh
    end
    
    properties (Access = private)
        nx1
        nx2
        coord
        perCoord
        epsilon
        nCells
    end
    
    methods (Access = public)
        
        function obj = Dehomogenizing()
            obj.init();
            obj.createCoord(); 
            obj.createBackgroundMesh();
            obj.createBoundaryMesh();                       
            obj.createUnfittedMesh();
            obj.createGrid(); 
            obj.createDistortedGrid();
            figure()
            hold on;                     
            obj.uMesh.plotStructureInColor([0.5 0.5 0.5]);
            obj.gMesh.plotStructureInColor('red');
            obj.dMesh.plotStructureInColor('black');
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.nx1    = 552;
            obj.nx2    = 552;
            obj.nCells = 15;
        end
        
        function createCoord(obj)
            x1 = linspace(0,1,obj.nx1);
            x2 = linspace(0,1,obj.nx2);
            x1T = repmat(x1,obj.nx2,1);
            x2T = repmat(x2',1,obj.nx1);
            obj.coord = [x1T(:),x2T(:)];
        end
        
        function createBackgroundMesh(obj)
            xy = obj.coord;                        
            x1 = xy(:,1);
            x2 = xy(:,2);
            connec   = delaunay(x1,x2);
            s.coord  = [x1,x2];
            s.connec = connec;
            obj.backgroundMesh = Mesh(s);  
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
            s.mMin  = 0.58;
            s.mMax  = 0.28;
            s.qMin  = 36;
            s.qMax  = 36;
            s.thetaMin = 0;
            s.thetaMax = -1*pi/16;
            ls = obj.createSuperEllipseLevelSet(s);
        end
        
        function ls = createSuperEllipseLevelSet(obj,s)
            s.coord = obj.coord;
            superEllipse = SuperEllipseDistributionExample(s);
            superEllipse.computeParameters();
            s.type  = 'periodicAndOriented';            
            s.coord = obj.coord;
            s.ndim  = 2;            
            s.angle = superEllipse.theta;  
            s.nCells = obj.nCells;
            s.cellLevelSetParams = obj.createLevelSetCellParams(superEllipse);
            levelSet = LevelSetCreator.create(s);            
            ls = levelSet.getValue();                        
        end     
        
        function createGrid(obj)
            ls = obj.createGridLevelSet();                                    
            s.boundaryMesh   = obj.boundaryMesh;
            s.backgroundMesh = obj.backgroundMesh;
            obj.gMesh = UnfittedMesh(s);
            obj.gMesh.compute(ls);
        end 
        
        function ls = createGridLevelSet(obj)  
            s.coord = obj.coord;
            s.mMin  = 0.95;
            s.mMax  = 0.95;
            s.qMin  = 36;
            s.qMax  = 36;
            s.thetaMin = 0;
            s.thetaMax = 0;
            ls = obj.createSuperEllipseLevelSet(s);
        end                
        
        function createDistortedGrid(obj)
            ls = obj.createDistortedGridLevelSet();                                    
            s.boundaryMesh   = obj.boundaryMesh;
            s.backgroundMesh = obj.backgroundMesh;
            obj.dMesh = UnfittedMesh(s);
            obj.dMesh.compute(ls);
        end    
        
        function ls = createDistortedGridLevelSet(obj)  
            s.coord = obj.coord;
            s.mMin  = 0.95;
            s.mMax  = 0.95;
            s.qMin  = 36;
            s.qMax  = 36;
            s.thetaMin = 0;
            s.thetaMax = -1*pi/16;
            ls = obj.createSuperEllipseLevelSet(s);
        end                           
        
    end
    
    methods (Access = private, Static)
        
        function s = createLevelSetCellParams(distr)
           s.type   = 'smoothRectangle';            
           q        = distr.q;             
           s.widthH = distr.m1;
           s.widthV = distr.m2;
           s.pnorm  = q;
           s.ndim   = 2;
        end   
                   
    end
    
end