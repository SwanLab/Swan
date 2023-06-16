classdef ExampleShapeDeriv < handle
    
    properties (Access = private)
       tolerance
       mesh 
       quadrature
       projector
       cost       
    end
    
    methods (Access = public)
        
        function obj = ExampleShapeDeriv()
            obj.init()
            obj.createInitialMesh();
            obj.createQuadrature();
            obj.createProjector();
            obj.compute()
        end
        
    end
    
    methods (Access = private)

        function init(obj)
            obj.tolerance = 1e-10;
        end
        
        function m = createInitialMesh(obj)
            file = 'CantileverBeam_Triangle_Linear';
            a.fileName = file;
            s = FemDataContainer(a);
            m = s.mesh;
            obj.mesh = m;
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

        function compute(obj)
            obj.cost(1) = obj.computeCost();
            tau = 0.02;
            hasNotConverged = true;
            iter = 2;
            while hasNotConverged 
                f  = obj.createFunction() ;
                df = obj.computeDerivative();                
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
        
        function remesh(obj)
            s.coord  = obj.mesh.coord;
            s.connec = delaunay(obj.mesh.coord);
            m = Mesh(s);
            obj.mesh = m;
        end
          
        
        function f = createFunction(obj)
            %%elipse%%
          %  s.fHandle = @(x) (x(1,:,:)-1).^2/1.5 + (x(2,:,:)-0.5).^2/0.5 -1;
            %%cor%%
            %  s.fHandle = @(x) ( (x(1,:,:)-1).^2 + (x(2,:,:)-0.3).^2 - 0.5 ).^3 - (x(1,:,:)-1).^2.*(x(2,:,:)-0.3).^3 ;
            %%flor%%
             s.fHandle = @(x) ((x(1,:,:)-1).^2 + (x(2,:,:)-0.5).^2 ).^3 - 4*(x(1,:,:)-1).^2.*(x(2,:,:)-0.5).^2 ;
            
            %%PDE%%
            % a = PDEShapeDerivative() ;
            % u = a.computeTemperature() ;
            % u_t = a.createTemperatureTarget() ;
            % s.fHandle = @(x) u - u_t ;
            
            
            s.ndimf   = 1;
            s.mesh    = obj.mesh;
            f = AnalyticalFunction(s);
        end
        
        function dfC = computeDerivative(obj)
            %%elipse%%
            %df{1} = @(x)  2*(x(1,:,:)-1)/1.5;
            %df{2} = @(x)  2*(x(2,:,:)-0.5)/0.5;
            %%cor%%
            df{1} = @(x)  6*( (x(1,:,:)-1).^2 + (x(2,:,:)-0.3).^2 - 0.5 ).^2.*(x(1,:,:)-1) - 2*(x(1,:,:)-1).*(x(2,:,:)-0.3).^3;
            df{2} = @(x)  6*( (x(1,:,:)-1).^2 + (x(2,:,:)-0.3).^2 - 0.5).^2.*(x(2,:,:)-0.3) - 3*(x(1,:,:)-1).^2.*(x(2,:,:)-0.3).^2;                        
            %%flor%%
            df{1} = @(x) 6*( (x(1,:,:)-1).^2 + (x(2,:,:)-0.5).^2 ).^2.*(x(1,:,:)-1) - 8*(x(1,:,:)-1).*(x(2,:,:)-0.5).^2 ;
            df{2} = @(x) 6*( (x(1,:,:)-1).^2 + (x(2,:,:)-0.5).^2 ).^2.*(x(2,:,:)-0.5) - 8*(x(2,:,:)-0.5).*(x(1,:,:)-1).^2 ;
            for i = 1:length(df)
                s.fHandle = df{i};
                s.ndimf   = 1;
                s.mesh    = obj.mesh;
                dfC{i} = AnalyticalFunction(s);
            end
        end
        
   
        function c = computeCost(obj)
            f    = obj.createFunction();
            fG   = f.evaluate(obj.quadrature.posgp);
            dVG  = obj.mesh.computeDvolume(obj.quadrature);
            fG   = squeeze(fG);
            c    = fG.*dVG;
            c = sum(c(:));
        end

        function plotCost(obj)
            figure(2)
            plot(obj.cost)
        end

        function addLegendToFigures(obj)
            figure(2)
            set(groot,'defaultAxesTickLabelInterpreter','latex');
            set(groot,'defaulttextinterpreter','latex');
            set(groot,'defaultLegendInterpreter','latex');
            
            fig1 = figure(2);
            hold on;
            title("\textbf{Cost vs number of iterations}");
            plot(c, 'b', 'LineWidth', 1);
            xlabel("$Number\ of\ iterations$");
            ylabel("$Cost\  of\  the\  function$ ");
            grid on;
            grid minor;
            box on;
            set(gcf, 'units', 'points', 'position', [50,50,676/2,420/2]);
            hold off;
            
            set(fig1, 'units', 'points');
            pos = get(fig1, 'position');
            set(fig1, 'PaperPositionMode', 'Auto', 'PaperUnits', 'points', ...
                'PaperSize',[pos(3), pos(4)]);
            
        end        
    end
    
end

