classdef PDEShapeDerivative < handle
    
    properties (Access = private)
        mesh
        temperatureTarget
        temperature
        adjoint

        quadrature
        tolerance
        cost
        projector

    end
    
    methods (Access = public)
        
        function obj = PDEShapeDerivative()
            addpath(genpath(fileparts(mfilename('fullpath'))))
            obj.init();  
            obj.createMesh();                        
            obj.createTemperatureTarget();
            obj.computeTemperature();
            obj.computeAdjoint();
            obj.createFunction() ;
            obj.createQuadrature();
            obj.createProjector();
            obj.compute()
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
            fem.print('asa') 
             obj.temperature = fem.temperature;
             % obj.temperature.plot()
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
            % obj.adjoint.plot()
        end

        function f = createRHSFunction(obj)
            s.f1 = obj.temperature;
            s.f2 = obj.temperatureTarget;
            s.operation = @(f1,f2) 2*f1 - f2;
            s.ndimf = 1 ;
            f = ComposedFunction(s);
        end     

        function t = createTemperatureTarget(obj)
            s.fHandle = @(x) 0.6^2 - (x(1,:,:)-0.5).^2 + (x(2,:,:)-0.5).^2 - 0.5;
            s.ndimf   = 1;
            s.mesh    = obj.mesh;
            t = AnalyticalFunction(s);  
            obj.temperatureTarget = t;
        end

        function f = createFunction(obj)
            %%PDE%%
            u = obj.temperature ;
            u_t = obj.temperatureTarget ;
            s.fHandle = @(x) (u - u_t).^2 ;
                        
            s.ndimf   = 1;
            s.mesh    = obj.mesh;
            f = AnalyticalFunction(s);
        end

        %%%%%%
        % function dfC = computeDerivative(obj)
        %     df{1} = @(x) 6*( (x(1,:,:)-1).^2 + (x(2,:,:)-0.5).^2 ).^2.*(x(1,:,:)-1) - 8*(x(1,:,:)-1).*(x(2,:,:)-0.5).^2 ;
        %     df{2} = @(x) 6*( (x(1,:,:)-1).^2 + (x(2,:,:)-0.5).^2 ).^2.*(x(2,:,:)-0.5) - 8*(x(2,:,:)-0.5).*(x(1,:,:)-1).^2 ;
        %     for i = 1:length(df)
        %         s.fHandle = df{i};
        %         s.ndimf   = 1;
        %         s.mesh    = obj.mesh;
        %         dfC{i} = AnalyticalFunction(s);
        %     end
        % end

        function c = computeCost(obj)
            f    = obj.createFunction();
            fG   = f.evaluate(obj.quadrature.posgp);
            dVG  = obj.mesh.computeDvolume(obj.quadrature);
            fG   = squeeze(fG);
            c    = fG.*dVG;
            c = sum(c(:));
        end

        function compute(obj)
            obj.cost(1) = obj.computeCost();
            tau = 0.02;
            hasNotConverged = true;
            iter = 2;
            while hasNotConverged 
                f  = obj.temperature ;
                df = obj.adjoint;                
                g = obj.projector.project(obj.mesh,f,df);
                obj.updateMesh(g,tau);
                cTrial =  obj.computeCost();                
                if cTrial < obj.cost(iter-1)
                    obj.cost(iter) = cTrial;
                    tau = 1.1*tau;
                    incC = abs(obj.cost(iter) - obj.cost(iter-1))/abs(obj.cost(iter));
                    hasNotConverged = incC > obj.tolerance;  
                    iter   = iter + 1 ;
                    obj.plotMesh();
                    obj.plotCost();
                   % obj.remesh();
                else 
                    tau  = tau/2;
                    if tau < obj.tolerance
                      hasNotConverged = false;
                    else
                      hasNotConverged = true;
                    end
                end
            end
          end

        function updateMesh(obj,g,tau)
            coord = obj.mesh.coord;
            nCoord(:,1) = coord(:,1) - g.fValues(:,1)*tau ;
            nCoord(:,2) = coord(:,2) - g.fValues(:,2)*tau ;
            s.connec = obj.mesh.connec;
            s.coord  = nCoord;
            m = Mesh(s);
            obj.mesh =  m;
        end

        function q = createQuadrature(obj)
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature('CUBIC');
            obj.quadrature = q;
        end     

        function createProjector(obj)
            p = ShapeDerivProjector();
            obj.projector = p;
        end
        
        function plotMesh(obj)
            figure(1)
            clf
            obj.mesh.plot()
            %newMesh = obj.remesh(newMesh);
            % axis([0 2 -0.3 1.4])
             axis([-0.25 2.25 -0.25 1.25])
            %axis([0 2 -0.5 1.3])
        end

        function plotCost(obj)
            figure(2)
            plot(obj.cost)
        end
        
    end
    
end