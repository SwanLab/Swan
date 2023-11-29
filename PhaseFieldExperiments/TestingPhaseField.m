classdef TestingPhaseField < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
        mesh
        boundaryConditions
    end
    
    methods (Access = public)
        
        function obj = TestingPhaseField()
            obj.init()
            obj.createMeshes();
            obj.solveProblem();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
           close all 
        end
        
        function solveProblem(obj)
            s.mesh = obj.mesh;
            p = PhaseFieldComputer(s);
            p.compute()            
        end
        
        function createMeshes(obj)
         %   obj.createOneElementMesh();
           %obj.createTwoElementMesh();
         %  obj.createArbitraryElementMesh();
            obj.createOpenHoleMesh();
        end       
        
        function createOneElementMesh(obj)        
            sM.coord = [0,0;
                1,0;
                1,1;
                0,1];
            sM.connec = [1 2 3 4];
            
            m = Mesh(sM);
            m.plot();
            obj.mesh = m;
        end        
        
        function createTwoElementMesh(obj)
            sM.coord = [-1,-1;
                1,-1;
                1,1;
                -1,1;
                -1 3;
                1 3];
            sM.connec = [1 2 3 4
                4 3 6 5];
            
            m = Mesh(sM);
            m.plot();           
            obj.mesh = m;
        end
        
       function createArbitraryElementMesh(obj)
            m = UnitQuadMesh(10, 10);
            obj.mesh = m;
        end        
        
        function createOpenHoleMesh(obj)
            %file = 'test2d_triangle';
            %a.fileName = file;
            % s = FemDataContainer(a);
            
            % Generate coordinates
            x1 = linspace(0,1,20);
            x2 = linspace(1,2,20);
            % Create the grid
            [xv,yv] = meshgrid(x1,x2);
            % Triangulate the mesh to obtain coordinates and connectiviti
            % x1 = linspace(0,1,20);
            x2 = linspace(0,1,20);
            [xv,yv] = meshgrid(x1,x2);
            [F,V] = mesh2tri(xv,yv,zeros(size(xv)),'x');
            sBg.coord  = V(:,1:2);
            sBg.connec = F;
            bgMesh = Mesh(sBg);
            bdMesh  = bgMesh.createBoundaryMesh();
            
            sLS.type       = 'circleInclusion';
            sLS.mesh       = bgMesh;
            sLS.ndim       = 2;
            sLS.fracRadius = 0.4;
            sLS.coord      = bgMesh.coord;
            ls = LevelSetCreator.create(sLS);
            levelSet = ls.getValue();
            
            sUm.backgroundMesh = bgMesh;
            sUm.boundaryMesh   = bdMesh;
            uMesh = UnfittedMesh(sUm);
            uMesh.compute(levelSet);
            
            obj.mesh = uMesh.createInnerMesh();
            obj.mesh = obj.mesh.computeCanonicalMesh();
            %obj.mesh = bgMesh;
            
            figure
            uMesh.plot
            [F,V] = mesh2tri(xv,yv,zeros(size(xv)),'x');
            
            s.coord = V(:,1:2);
            s.connec = F;
            m = Mesh(s);
            obj.mesh = m;
            
            x1 = linspace(0,1,30);
            x2 = linspace(0,1,30);
            [xv,yv] = meshgrid(x1,x2);
            [F,V] = mesh2tri(xv,yv,zeros(size(xv)),'x');
            sBg.coord  = V(:,1:2);
            sBg.connec = F;
            bgMesh = Mesh(sBg);
            bdMesh  = bgMesh.createBoundaryMesh();
            
            sLS.type       = 'circleInclusion';
            sLS.mesh       = bgMesh;
            sLS.ndim       = 2;
            sLS.fracRadius = 0.4;
            sLS.coord      = bgMesh.coord;
            ls = LevelSetCreator.create(sLS);
            levelSet = ls.getValue();
            
            sUm.backgroundMesh = bgMesh;
            sUm.boundaryMesh   = bdMesh;
            uMesh = UnfittedMesh(sUm);
            uMesh.compute(levelSet);
            
            %obj.mesh = uMesh.createInnerMesh();
            obj.mesh = uMesh.createInnerMeshGoodConditioning();
            %obj.mesh = obj.mesh.computeCanonicalMesh();
            obj.mesh.plot;

            
        end
        

        
    end
    
end