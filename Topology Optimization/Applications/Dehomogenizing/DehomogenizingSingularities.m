classdef DehomogenizingSingularities < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
        mesh
        orientation
    end
    
    methods (Access = public)
        
        function obj = DehomogenizingSingularities()
            obj.createMesh();
            obj.createOrientation();
            obj.computeSingularities();
        end
        
    end
    
    methods (Access = private)
        
        
        function createMesh(obj)
            obj.createBenchmarkMesh();            
        end
        
        function createOrientation(obj)
            obj.createBenchmarkOrientation();
        end
        
        function computeSingularities(obj)
            s.mesh        = obj.mesh;
            s.orientation = obj.orientation;
            sF = SingularitiesFinder(s);
            isS = sF.computeSingularElements();
            sF.plot();
        end
        
        function createBenchmarkMesh(obj)
           h = 0.03;
           xmin = 0.5;
           xmax = 2;
           ymin = 0.25;
           ymax = 1.75;
           xv = xmin:h:xmax;
           yv = ymin:h:ymax;
           [X,Y] = meshgrid(xv,yv);
           s.coord(:,1) = X(:);
           s.coord(:,2) = Y(:);
           s.connec = delaunay(s.coord); 
%             [F,V] = mesh2tri(X,Y,zeros(size(X)),'x');
%             s.coord  = V(:,1:2);
%               s.connec = F;

            m = Mesh(s);            
            %m.plot;     
            obj.mesh = m;

%             s.coord = [0 0; 0 1; 1 1; 1 0];
%             s.connec = [1 2 4; 2 3 4];
%             m = Mesh(s);
%             obj.mesh = m;
        end        
        
        function createBenchmarkOrientation(obj)
            s1 = 0.3;
            s2 = -0.8;
            x1 = obj.mesh.coord(:,1);
            x2 = obj.mesh.coord(:,2);
            v(:,1) = cos(pi*(x1 + s1*x2));
            v(:,2) = cos(pi*(x2 + s2*x1));
            nv  = sqrt(v(:,1).^2 + v(:,2).^2);
            ex  = [1 0];
            nex = sqrt(ex(:,1).^2 + ex(:,2).^2);
            vex = v(:,1)*ex(:,1) + v(:,2)*ex(:,2);
            beta = acos(vex./(nv.*nex));
            alpha = beta/2;
            obj.orientation(:,1) = cos(alpha);
            obj.orientation(:,2) = sin(alpha);
        end        
        
        
        
    end
    
end