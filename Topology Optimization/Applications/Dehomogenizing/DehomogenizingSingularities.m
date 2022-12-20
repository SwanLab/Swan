classdef DehomogenizingSingularities < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
        mesh
        orientation
        backgroundMesh
    end
    
    methods (Access = public)
        
        function obj = DehomogenizingSingularities()
            obj.createMesh();
            obj.createBackgroundMesh();            
            obj.createOrientation();
            obj.dehomogenize();
        end
        
    end
    
    methods (Access = private)
        
        
        function createMesh(obj)
            obj.createBenchmarkMesh();            
        end
        
        function createOrientation(obj)
            m = obj.backgroundMesh;%obj.mesh
            s1 = 0.32;
            s2 = -0.8;
            x1 = m.coord(:,1);
            x2 = m.coord(:,2);
            v(:,1) = cos(pi*(x1 + s1*x2));
            v(:,2) = cos(pi*(x2 + s2*x1));
            beta = atan2(v(:,2),v(:,1));
            alpha = beta/2;
            obj.orientation(:,1) = cos(alpha);
            obj.orientation(:,2) = sin(alpha);
            
            
           figure()
            a = obj.orientation;
            x = obj.mesh.coord(:,1);
            y = obj.mesh.coord(:,2);
            ax = a(:,1);
            ay = a(:,2);
            q = quiver(x1,x2,ax,ay);
            q.ShowArrowHead = 'off';
        end
        
        function createBenchmarkMesh(obj)
           h = 0.0221;
           xmin = 0.50;
           xmax = 2.0;
           ymin = 0.25;
           ymax = 1.75;
           xv = xmin:h:xmax;
           yv = ymin:h:ymax;
           [X,Y] = meshgrid(xv,yv);
           s.coord(:,1) = X(:);
           s.coord(:,2) = Y(:);
           s.connec = delaunay(s.coord); 
%              [F,V] = mesh2tri(X,Y,zeros(size(X)),'x');
%              s.coord  = V(:,1:2);
%                s.connec = F;
            m = Mesh(s);            
            obj.mesh = m;
        end        
               

        function dehomogenize(obj)
            s.backgroundMesh     = obj.backgroundMesh;
            s.nCells             = 45;
            s.cellLevelSetParams = obj.createLevelSetCellParams();
            s.mesh               = obj.backgroundMesh;%obj.mesh;
            s.theta              = atan2(obj.orientation(:,2),obj.orientation(:,1));
            d = Dehomogenizer(s);
            d.compute();
            d.plot();
        end


         function s = createLevelSetCellParams(obj)        
            s.type   = 'rectangleInclusion';
            s.widthH = 0.87*ones(size(obj.backgroundMesh.coord,1),1);
            s.widthV = 0.87*ones(size(obj.backgroundMesh.coord,1),1);        
            s.ndim   = 2;
         end          

         function createBackgroundMesh(obj)
             FV.vertices = [obj.mesh.coord,zeros(size(obj.mesh.coord,1),1)];
             FV.faces    = obj.mesh.connec;
             FV2 = FV;
             FV2 = refinepatch(FV2);
             FV2 = refinepatch(FV2);
          %   FV2 = refinepatch(FV2);
           %  FV2 = refinepatch(FV2);             
             s.coord = FV2.vertices(:,1:2);
             s.connec = FV2.faces;
             m = Mesh(s);
             obj.backgroundMesh = m;
         end                
        
        
    end
    
end