classdef Tutorial00_2_Remeshing < handle
    

    properties (Access = private)
        mesh
        analyticalFunction
    end

    properties (Access = private)
        nLevels        
    end    
    
    methods (Access = public)
        
        function obj = Tutorial00_2_Remeshing()
            obj.init();
            obj.createMesh();
            obj.createAnalyticalFunction();
            obj.remesh();
        end

    end
    
    methods (Access = private)

        function init(obj)
            obj.nLevels = 3;
        end
                
        function createMesh(obj)
            xmin = 0; 
            xmax = 1;
            ymin = 0;
            ymax = 1;
            h = 0.1;
            xv = xmin:h:xmax;
            yv = ymin:h:ymax;
            [X,Y] = meshgrid(xv,yv);
            s.coord(:,1) = X(:);
            s.coord(:,2) = Y(:);
            [F,V] = mesh2tri(X,Y,zeros(size(X)),'f');
            s.coord  = V(:,1:2);
            s.connec = F;
            m = Mesh.create(s);
            obj.mesh = m;
            obj.mesh.plot();             
        end

        function f = createAnalyticalFunction(obj)
            s.fHandle = @(x) obj.circle(x);
            s.ndimf   = 2;
            s.mesh    = obj.mesh;
            f       = AnalyticalFunction(s);
            obj.analyticalFunction = f;
        end

        function remesh(obj)
            f = obj.analyticalFunction();
            fC = f.project('P1');
            fC.plot();
            fD = f.project('P1D');
            fD.plot()
            m = obj.mesh;
            for iLevel = 1:obj.nLevels
                m  = m.remesh();
                fC = fC.refine(m);
                fC.plot();
                fD = fD.refine(m);
                fD.plot();
            end
        end

    end

    methods (Static, Access = private)
        function f = circle(x)
            x1 = x(1,:,:);
            x2 = x(2,:,:);
            f1 = 1-heaviside((x1-0.5).^2+(x2-0.5).^2-0.3.^2);
            f2 = 1-heaviside((x1-0.1).^2+(x2-0.1).^2-0.7.^2);
           f = [f1;f2];
        end
        
    end
    
end