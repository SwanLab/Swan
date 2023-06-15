classdef PDEShapeDerivative < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        mesh
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = PDEShapeDerivative()
            obj.init();
            obj.createMesh();
            obj.createFEM();
           % obj.createAdjoint();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            
        end
        
        function createMesh(obj)
             x1 = linspace(0,1,10);
             x2 = linspace(0,1,10);             
             [xv,yv] = meshgrid(x1,x2); 
             [F,V] = mesh2tri(xv,yv,zeros(size(xv)),'x');
             s.coord  = V(:,1:2);
             s.connec = F;
             obj.mesh = Mesh(s);         
        end

        function bc = createDirichletData(obj)
           bc.dirichlet(:,1) = obj.createBoundaryNodes();
           bc.dirichlet(:,2) = 0;
           bc.pointload = [];
        end

        function nodesB = createBoundaryNodes(obj)
           boundaryMesh = obj.mesh.createBoundaryMesh();
           for iFace = 1:numel(boundaryMesh)
                bM = boundaryMesh{iFace};
                nodesI(:,iFace) = bM.nodesInBoxFaces;
           end
           nodesB = any(nodesI,2);        
        end        
        
        function createFEM(obj)
            s.bc    = obj.createDirichletData(); 
            s.scale = 'MACRO';
            s.type  = 'THERMAL';
            s.mesh  = obj.mesh;
            s.dim   = '2D';
            fem = FEM.create(s);
            fem.solve();            
        end

    end
    
end