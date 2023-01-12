classdef DehomogenizingSingularitiesTest < handle
        
    properties (Access = private)
        mesh
        orientation
        backgroundMesh
        levelSet
    end

    properties (Access = private)
        testName
        meshSize
        singularitiesData
        xmin
        xmax
        ymin
        ymax
        widthH
        widthW
        nCells
    end
    
    methods (Access = public)
        
        function obj = DehomogenizingSingularitiesTest(cParams)
            obj.init(cParams);
            obj.createMesh();
            obj.createBackgroundMesh();            
            obj.createOrientation();
            obj.dehomogenize();
        end

        function passed = hasPassed(obj)
            d = load(obj.testName);
            ls = obj.levelSet;
            itIs = isequaln(ls,d.levelSet);
            passed = itIs;
        end

    end
    
    methods (Access = private)

        function init(obj,cParams)
            obj.testName = cParams.testName;
            %obj.meshSize = 0.00521;
            obj.meshSize = 0.09;%0.0221;%0.0521 %0.0221;0.0921
            obj.nCells   = linspace(60,62,40);%45;   %45        
            obj.xmin = 0.5;
            obj.xmax = 2.0;
            obj.ymin = 0.25;
            obj.ymax = 1.75;
            obj.singularitiesData = [0.32,-0.8];            
            obj.widthH = 0.87;
            obj.widthW = 0.87;
        end        
        
        function createMesh(obj)
           h = obj.meshSize;
           xv = obj.xmin:h:obj.xmax;
           yv = obj.ymin:h:obj.ymax;
           [X,Y] = meshgrid(xv,yv);
           s.coord(:,1) = X(:);
           s.coord(:,2) = Y(:);
           s.connec = delaunay(s.coord);
           m = Mesh(s);
           obj.mesh = m;
        end

         function createBackgroundMesh(obj)
             FV.vertices = [obj.mesh.coord,zeros(size(obj.mesh.coord,1),1)];
             FV.faces    = obj.mesh.connec;
             FV2 = FV;
                FV2 = refinepatch(FV2);
             s.coord = FV2.vertices(:,1:2);
             s.connec = FV2.faces;
             m = Mesh(s);
             obj.backgroundMesh = m;
         end          
        
        function createOrientation(obj)
            m = obj.backgroundMesh;
            s1 = obj.singularitiesData(:,1);
            s2 = obj.singularitiesData(:,2);
            x1 = m.coord(:,1);
            x2 = m.coord(:,2);
            v(:,1) = cos(pi*(x1 + s1*x2));
            v(:,2) = cos(pi*(x2 + s2*x1));
            beta = atan2(v(:,2),v(:,1));
            alpha = beta/2;
            obj.orientation(:,1) = cos(alpha);
            obj.orientation(:,2) = sin(alpha);
        end

         function s = createLevelSetCellParams(obj)        
            I        = ones(size(obj.backgroundMesh.coord,1),1);
            s.type   = 'rectangleInclusion';
            s.widthH = obj.widthH*I;
            s.widthV = obj.widthW*I;
            s.ndim   = 2;
         end       

        function ls = dehomogenize(obj)
            s.backgroundMesh     = obj.backgroundMesh;
            s.nCells             = obj.nCells;
            s.cellLevelSetParams = obj.createLevelSetCellParams();
            s.mesh               = obj.backgroundMesh;
            s.theta              = atan2(obj.orientation(:,2),obj.orientation(:,1));
            d = Dehomogenizer(s);
            ls = d.compute();
            obj.levelSet = ls;
            d.plot();
        end         
        
    end
    
end