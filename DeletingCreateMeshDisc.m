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
            obj.createMesh();
            obj.mesh.plot(); 
            m = obj.mesh.createDiscontinuousMesh();
            fC = obj.createP1Function();
            fD = fC.project('P1D');
            fD.plot()
        end
        
    end
    
    methods (Access = private)
        
        function createMesh(obj)
            xmin = 0; 
            xmax = 1;
            ymin = 0;
            ymax = 1;
            h = 1;
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

        function f = createP1Function(obj)
            %s.fHandle = @(x) sin(10*x(1,:,:));
            s.fHandle = @(x) x;
            s.ndimf   = 1;
            s.mesh    = obj.mesh;            
            f = AnalyticalFunction(s);
            f = f.project('P1');
            f.plot()
        end
        
    end
    
end