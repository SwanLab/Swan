classdef DehomogenizingRadialExample < handle
    
    properties (Access = private)
        backgroundMesh
        theta        
        superEllipse
        cellLevelSetParams
    end
    
    properties (Access = private)
        nx1
        nx2
        nCells
    end
    
    methods (Access = public)
        
        function obj = DehomogenizingRadialExample()
            obj.init();
            obj.createBackgroundMesh();  
            obj.createOrientation();            
            obj.createSuperEllipseParams();
            obj.createLevelSetCellParams();
            obj.dehomogenize();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.nx1    = 125*2;
            obj.nx2    = 125;
            obj.nCells = 16;
        end
               
        function createBackgroundMesh(obj)
            
            x1 = linspace(-1,1,obj.nx1);
            x2 = linspace(0,1,obj.nx2);
           
            x1T = repmat(x1,obj.nx2,1);
            x2T = repmat(x2',1,obj.nx1);     
            
%             xy = obj.coord;                        
%             x1 = xy(:,1);
%             x2 = xy(:,2);
%             s.connec   = delaunay(x1,x2);
%             s.coord    = obj.coord;
%             obj.backgroundMesh = Mesh(s);
%             
                       
                             
%             x1min = min(x1);
%             x1max = max(x1);
%             x2min = min(x2);
%             x2max = max(x2);
%             [coordinates, nodes,nel,nnode] = MeshRectangular(x1max-x1min,x2max-x2min,obj.nx1,obj.nx2);
%             s.coord(:,1) = coordinates(:,1)+x1min;
%             s.coord(:,2) = coordinates(:,2)+x2min;
%             s.connec = nodes;
%             obj.backgroundMesh = Mesh(s);              
%             obj.backgroundMesh.plot();
            
                                
             [xv,yv] = meshgrid(x1,x2); 
             [F,V] = mesh2tri(xv,yv,zeros(size(xv)),'x');             
             s.coord  = V(:,1:2);
             s.connec = F;
             obj.backgroundMesh = Mesh(s);  
             obj.backgroundMesh.plot()
%             obj.coord = s.coord;
            
        end
        
        function createOrientation(obj)
            x2 = obj.backgroundMesh.coord(:,2);
            x1 = obj.backgroundMesh.coord(:,1);
            obj.theta = atan2(x2,x1);
        end        

        function createSuperEllipseParams(obj)
           s.coord = obj.backgroundMesh.coord;
           s.mMin  = 0.4;
           s.mMax  = 0.99;
           s.qMin  = 32;
           s.qMax  = 32;            
           sE = SuperEllipseDistributionExample(s);
           sE.computeParameters();            
           obj.superEllipse = sE;
        end
        
        function createLevelSetCellParams(obj)
           s.type   = 'smoothRectangle';            
           s.widthH = obj.superEllipse.m1;
           s.widthV = obj.superEllipse.m2;
           s.pnorm  = obj.superEllipse.q;
           s.ndim   = 2;
           obj.cellLevelSetParams = s;
        end   
        
        function dehomogenize(obj)
            s.backgroundMesh     = obj.backgroundMesh;
            s.nCells             = obj.nCells;
            s.theta              = obj.theta;
            s.cellLevelSetParams = obj.cellLevelSetParams;
            s.mesh               = obj.backgroundMesh;
            d = Dehomogenizer(s);
            d.compute();
            d.plot();            
        end        
                   
    end
    
end