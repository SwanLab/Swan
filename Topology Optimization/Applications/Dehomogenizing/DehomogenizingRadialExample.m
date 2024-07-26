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
            for i = 20:50
                obj.nCells = i;
                obj.createBackgroundMesh();
                obj.createOrientation();
                obj.createSuperEllipseParams();
                obj.createLevelSetCellParams();
                obj.dehomogenize();
                exportgraphics(gcf,'testAnimated.gif','Append',true);
            end
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.nx1    = 55*2;%180
            obj.nx2    = 55;%180
            obj.nCells = 42;%32
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
             obj.backgroundMesh = Mesh.create(s);
        %     obj.backgroundMesh.plot()
%             obj.coord = s.coord;
            
        end
        
        function createOrientation(obj)
            x2 = obj.backgroundMesh.coord(:,2);
            x1 = obj.backgroundMesh.coord(:,1);
            obj.theta = atan((x2+0.4)./x1);            
            isLeft = x1 < 0;
            %obj.theta(isLeft) = obj.theta(isLeft) + 180;
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
           s.type   = 'rectangleInclusion';%'smoothRectangle';
           s.widthH = 0.95*ones(size(obj.superEllipse.m1));
           s.widthV = 0.95*ones(size(obj.superEllipse.m2));
           s.pnorm  = obj.superEllipse.q;
           s.ndim   = 2;
           obj.cellLevelSetParams = s;
        end
        
        function dehomogenize(obj)

            alphaM(:,2) = cos(obj.theta);
            alphaM(:,1) = sin(obj.theta);
            

            s.fValues = alphaM;
            s.mesh    = obj.backgroundMesh;
            s.order   = 'P1';
            a1{1} = LagrangianFunction(s);
            %a1{1} = a0{1}.project('P1');


            s.fValues(:,1) = -alphaM(:,2);
            s.fValues(:,2) = alphaM(:,1);
            s.mesh    = obj.backgroundMesh;
            s.order   = 'P1';
            a1{2} = LagrangianFunction(s);
            %a1{2} = a0{2}.project('P1');


            s.backgroundMesh     = obj.backgroundMesh;
            s.nCells             = obj.nCells;
            s.theta              = a1;
            s.cellLevelSetParams = obj.cellLevelSetParams;
            s.mesh               = obj.backgroundMesh;
            d = Dehomogenizer(s);
            d.compute();
            d.plot();
        end
 
    end
    
end