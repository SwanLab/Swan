classdef Dehomogenizer < handle
    
    properties (Access = public)
       
    end
    
    properties (Access = private)
        boundaryMesh
        uMesh  
        remesher
        dilation
        phi
        epsilon
        levelSet
    end
    
    properties (Access = private)
        cellLevelSetParams
        backgroundMesh
        nCells
        theta
        mesh
    end
    
    methods (Access = public)
        
        function obj = Dehomogenizer(cParams)
            obj.init(cParams)            
        end

        function createRemesher(obj)
            m = obj.backgroundMesh;
            m = m.createDiscontinuousMesh();
            s.mesh = m;
            s.nLevels = 2;
            r  = Remesher(s);
            r.remesh();
            obj.remesher = r;
        end
        
        function ls = compute(obj)
            obj.computeDilation();
            obj.createMapping();
            obj.createRemesher();
            nCell = linspace(20,22,40);
            for i = 1:length(nCell)
                obj.nCells = nCell(i);
                obj.createEpsilon();
                obj.createLevelSet();                                    
                obj.plot();
                xmin = min(obj.backgroundMesh.coord(:,1));
                xmax = max(obj.backgroundMesh.coord(:,1));
                ymin = min(obj.backgroundMesh.coord(:,2));
                ymax = max(obj.backgroundMesh.coord(:,2));
                axis([xmin xmax ymin ymax])
                set(gca, 'Visible', 'off')                
                exportgraphics(gcf,'testAnimated2.gif','Append',true);
                close all
            end
            ls = obj.levelSet;            
        end
        
        function plot(obj)
            fineMesh = obj.remesher.fineMesh;
            m = fineMesh.createDiscontinuousMesh();
            obj.backgroundMesh = m;
            obj.createBoundaryMesh();                
            obj.createUnfittedMesh();                        
         %   obj.plotOrientation();
            obj.plotStructure();
           % obj.plotComponents();
       %    obj.backgroundMesh = m0;
        %   obj.createBoundaryMesh();                           
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.backgroundMesh     = cParams.backgroundMesh;
            obj.nCells             = cParams.nCells;
            obj.theta              = cParams.theta;
            obj.cellLevelSetParams = cParams.cellLevelSetParams;            
            obj.mesh               = cParams.mesh;
        end
        
        function createBoundaryMesh(obj)
            sB.backgroundMesh = obj.backgroundMesh;
            sB.dimension = 1:3;
            sB.type = 'FromReactangularBox';
            bMc = BoundaryMeshCreator.create(sB);
            obj.boundaryMesh  = bMc.create();            
        end
        
        function createUnfittedMesh(obj)
            ls = obj.levelSet;
            s.boundaryMesh   = obj.boundaryMesh;
            s.backgroundMesh = obj.backgroundMesh;
            obj.uMesh = UnfittedMesh(s);
            obj.uMesh.compute(ls);
        end               

        function computeDilation(obj)
            s.theta = obj.theta;
            s.mesh  = obj.mesh;
            dC = DilationFieldComputer(s);
            d  = dC.compute();
            %dC.plot();
            obj.dilation = d;
        end      

        function createMapping(obj)
            s.mesh     = obj.mesh;
            s.theta    = obj.theta;
            s.dilation = obj.dilation;
            c = ConformalMappingComputer(s);
            phiV = c.compute();
            % c.plot();
            obj.phi = phiV;
        end

        function createEpsilon(obj)
            L = obj.mesh.computeCharacteristicLength();
            obj.epsilon = L/obj.nCells;            
        end            
        
        function plotOrientation(obj)
            figure()
            x = obj.mesh.coord(:,1);
            y = obj.mesh.coord(:,2);
            t  = obj.theta;
            ct = cos(t(:,1));
            st = sin(t(:,1));
            quiver(x,y,ct,st)
        end
        
        function plotStructure(obj)
            figure()
            obj.uMesh.plotStructureInColor('black');
        end
        
        function plotComponents(obj)
            s.unfittedMesh = obj.uMesh;
            sp = UnfittedMeshSplitter(s);
            sp.split();                        
            sp.plot();
        end

        function createLevelSet(obj)
            s.coord  = obj.backgroundMesh.coord;            
            s.type   = 'periodicAndOriented';            
            s.backgroundMesh   = obj.backgroundMesh;
            s.mesh   = obj.mesh;
            s.remesher = obj.remesher;
            s.ndim   = 2;            
            s.phi = obj.phi;            
            s.dilation = obj.dilation;
            s.epsilon = obj.epsilon;
            s.cellLevelSetParams = obj.cellLevelSetParams;
            lSet = LevelSetCreator.create(s);            
            obj.levelSet = lSet.getValue();   
        end             
        
    end
    
end