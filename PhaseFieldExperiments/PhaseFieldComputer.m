classdef PhaseFieldComputer < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
        mesh
        boundaryConditions
        materialInterpolation
    end
    
    methods (Access = public)
        
        function obj = PhaseFieldComputer()
            obj.init()
            obj.computeMesh();
            obj.createBoundaryConditions();
            obj.createMaterialInterpolation();
            obj.computeFEM();
        end
        
       
        
        
    end
    
    
    
    methods (Access = private)
        
        function init(obj)
            
        end
        
        function computeMesh(obj)
            %file = 'test2d_triangle';
            %a.fileName = file;
            % s = FemDataContainer(a);
                       
            % Generate coordinates
            x1 = linspace(0,2,20);
            x2 = linspace(1,2,20);
            % Create the grid
            [xv,yv] = meshgrid(x1,x2);
            % Triangulate the mesh to obtain coordinates and connectivities
            [F,V] = mesh2tri(xv,yv,zeros(size(xv)),'x');
            
            s.coord = V(:,1:2);
            s.connec = F;
            m = Mesh(s);
            obj.mesh = m;
        end
        
      function createBoundaryConditions(obj)
            dirichletNodes = abs(obj.mesh.coord(:,1)-0) < 1e-12;
            rightSide  = max(obj.mesh.coord(:,1));
            isInRight = abs(obj.mesh.coord(:,1)-rightSide)< 1e-12;
            isInMiddleEdge = abs(obj.mesh.coord(:,2)-1.5) < 0.1;
            forceNodes = isInRight & isInMiddleEdge;
            nodes = 1:obj.mesh.nnodes;
            bc.dirichlet = nodes(dirichletNodes);
            bc.pointload(:,1) = nodes(forceNodes);
            bc.pointload(:,2) = 2;
            bc.pointload(:,3) = -1;
            obj.boundaryConditions = bc;
      end
        
      function createPhaseField()
          
      end
        
        function createMaterialInterpolation(obj)
            
            
            

            
           
            
            
        end
        
            
        function  createMaterial(obj,mesh,ngaus)
            int = obj.materialInterpolation;
            I = ones(mesh.nelem,ngaus);
            s.ptype = 'ELASTIC';
            s.pdim  = '2D';
            s.nelem = mesh.nelem;
            s.mesh  = mesh;
            s.kappa = int.kappa;%.9107*I;
            s.mu    = int.mu;%.3446*I;
            mat = Material.create(s);
            mat.compute(s);
            obj.material = mat;
        end
   
        
        function computeFEM(obj)
            s.mesh = obj.mesh;
            s.type = 'ELASTIC';
            s.scale = 'MACRO';
            s.material = obj.createMaterial(obj.mesh,1);
            s.dim = '2D';
            s.bc = obj.boundaryConditions;
            fem = FEM.create(s);
            fem.solve();
            
            figure(1)
            fem.uFun.plot()
            figure(2)
            fem.stressFun.plot()
            figure(3)
            fem.strainFun.plot()            
        end
        
        

        
    end
    
end
