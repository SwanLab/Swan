classdef PDEShapeDerivative < handle
    
    properties (Access = private)
        mesh
        temperature
        adjoint
    end
    
    methods (Access = public)
        
        function obj = PDEShapeDerivative()
            obj.init();            
            obj.createMesh();            
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

        function rhs = createVolumetricRHS(obj)
            s.fHandle = @(x,y) 4;
            s.ndimf   = 1;
            s.mesh    = obj.mesh;
            rhs = AnalyticalFunction(s);                         
        end
        
        function computeTemperature(obj)
            s.bc     = obj.createDirichletData(); 
            s.scale  = 'MACRO';
            s.type   = 'THERMAL';
            s.mesh   = obj.mesh;
            s.rhsVol = obj.createVolumetricRHS();
            fem = FEM.create(s);
            fem.solve(); 
            %fem.print('asa')
            obj.temperature = fem.temperature;
        end

        function computeAdjoint(obj)
            p1 = obj.computeGeneralAdjoint(obj.temperature);
            p2 = obj.computeGeneralAdjoint(obj.createTemperatureTarget());
            p = -2*(p1.fValues - p2.fValues);
            s.fValues = p;
            s.mesh    = obj.mesh;
            p = P1Function(s);
            obj.adjoint = p;
        end

        function p = computeGeneralAdjoint(obj,fun)
            s.bc     = obj.createDirichletData(); 
            s.scale  = 'MACRO';
            s.type   = 'THERMAL';
            s.mesh   = obj.mesh;
            s.rhsVol = fun;
            fem = FEM.create(s);
            fem.solve();  
            p = fem.temperature;
        end        

        function p = comptueAdjoint2(obj)
            s.bc     = obj.createDirichletData(); 
            s.scale  = 'MACRO';
            s.type   = 'THERMAL';
            s.mesh   = obj.mesh;
            s.rhsVol = obj.createTemperatureTarget();
            fem = FEM.create(s);
            fem.solve(); 
            p = fem.temperature;            
        end        

        function rhs = createTemperatureTarget(obj)
            s.fHandle = @(x) 0.6^2 - (x(:,1)-0.5).^2 + (x(:,2)-0.5).^2 - 0.5;
            s.ndimf   = 1;
            s.mesh    = obj.mesh;
            rhs = AnalyticalFunction(s);  
        end



    end
    
end