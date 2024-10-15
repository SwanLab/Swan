classdef DeletingCreateMeshDisc < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
        mesh
    end
    
    methods (Access = public)
        
        function obj = DeletingCreateMeshDisc()
            close all
            obj.createMesh();
            obj.mesh.plot(); 
      %      m = obj.mesh.createDiscontinuousMesh();
            f = obj.createAnalyticalFunction(obj.mesh);
            fC = f.project('P1');
            fC.plot();
            %d = Grad(fC);
            %d.plot(obj.mesh);
            fD = f.project('P1D');
            fD.plot()

            fD2 = fC.project('P1');

            
            m = obj.mesh;
            for i = 1:3
                m  = m.remesh();
                %fi  = obj.createAnalyticalFunction(m);
                %fCi = fi.project('P1'); fCi.plot
                fC = fC.refine(m); 
                fC.plot();    

                fD = fD.refine(fD,m);
                fD.plot();
            end
        
        end
        
    end
    
    methods (Access = private)
        
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
        end

        function f = createAnalyticalFunction(obj,m)

            s.fHandle = @(x) obj.circle(x);%x(1,:,:);
            s.ndimf   = 1;
            s.mesh    = m;
            f       = AnalyticalFunction(s);
            
            % %s.fHandle = @(x) sin(10*x(1,:,:));
            % s.fHandle = @(x) x.^(x-1);
            % s.ndimf   = 1;
            % s.mesh    = obj.mesh;            
            % f = AnalyticalFunction(s);


        end
    end

    methods (Static, Access = private)
        function f = circle(x)
            x1 = x(1,:,:);
            x2 = x(2,:,:);
            f1 = 1-heaviside((x1-0.5).^2+(x2-0.5).^2-0.3.^2);
         %   f2 = 1-heaviside((x1-0.1).^2+(x2-0.1).^2-0.7.^2);
          %  f = [f1;f2];
            f = f1;
        end
        
    end
    
end