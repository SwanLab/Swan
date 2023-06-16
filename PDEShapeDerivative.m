classdef PDEShapeDerivative < handle
    
    properties (Access = private)
        mesh
        temperatureTarget
        temperature
        adjoint
    end
    
    methods (Access = public)
        
        function obj = PDEShapeDerivative()
            obj.init();  
            obj.createMesh();                        
            obj.createTemperatureTarget();
            obj.computeTemperature();
            obj.computeAdjoint();
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
           bc.dirichlet(:,2) = 1;
           bc.dirichlet(:,3) = 0;
           bc.pointload = [];
        end

        function nodesB = createBoundaryNodes(obj)
           boundaryMesh = obj.mesh.createBoundaryMesh();
           for iFace = 1:numel(boundaryMesh)
                bM = boundaryMesh{iFace};
                nodesI(:,iFace) = bM.nodesInBoxFaces;
           end
           nodesB = any(nodesI,2);  
           nodesB = find(nodesB);
        end     

        function rhs = createFunctionRHS(obj)
            s.fHandle = @(x,y) 4*ones(size(x));
            s.ndimf   = 1;
            s.mesh    = obj.mesh;
            rhs = AnalyticalFunction(s);                         
        end
        
        function computeTemperature(obj)
            s.bc     = obj.createDirichletData(); 
            s.scale  = 'MACRO';
            s.type   = 'THERMAL';
            s.mesh   = obj.mesh;
            s.rhsFun = obj.createFunctionRHS();
            fem = FEM.create(s);
            fem.solve(); 
            %fem.print('asa')
            obj.temperature = fem.temperature;
        end

        function computeAdjoint(obj)
            s.bc     = obj.createDirichletData(); 
            s.scale  = 'MACRO';
            s.type   = 'THERMAL';
            s.mesh   = obj.mesh;
            s.rhsFun = obj.createRHSFunction();
            fem = FEM.create(s);
            fem.solve();  
            p = fem.temperature;
            obj.adjoint = p;
        end

        function f = createRHSFunction(obj)
            s.f1 = obj.temperature;
            s.f2 = obj.temperatureTarget;
            s.operation = @(f1,f2) 2*f1 - f2;
            f = ComposedFunction(s);
        end     

        function t = createTemperatureTarget(obj)
            s.fHandle = @(x) 0.6^2 - (x(1,:,:)-0.5).^2 + (x(2,:,:)-0.5).^2 - 0.5;
            s.ndimf   = 1;
            s.mesh    = obj.mesh;
            t = AnalyticalFunction(s);  
            obj.temperatureTarget = t;
        end



    end
    
end