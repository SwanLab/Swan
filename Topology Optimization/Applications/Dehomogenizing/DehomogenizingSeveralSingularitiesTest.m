classdef DehomogenizingSeveralSingularitiesTest < handle
        
    properties (Access = private)
        mesh
        orientation
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
        
        function obj = DehomogenizingSeveralSingularitiesTest(cParams)
            obj.init(cParams);
            obj.createMesh();
            obj.createOrientation();
            obj.dehomogenize();
        end

        function passed = hasPassed(obj)
            d = load(obj.testName);
            ls = obj.levelSet;
            errI = zeros(numel(ls),1);
            for iCell = 1:numel(ls)
                lSI = d.levelSet{iCell};
                errI(iCell) = norm(ls{iCell} - lSI)/norm(lSI);
            end
            itIs = sum(errI) < 1e-12;
            passed = itIs;
        end

    end
    
    methods (Access = private)

        function init(obj,cParams)
            obj.testName = cParams.testName;
            %obj.meshSize = 0.00521;
            obj.meshSize = 0.09;%0.0221;%0.09;%0.0221;%0.0521 %0.0221;0.0921
            obj.nCells   = [60 62];%linspace(60,62,40);%45;   %45        
            obj.xmin = 0.8;
            obj.xmax = 4.8;
            obj.ymin = 0.1;
            obj.ymax = 4.1;
            obj.singularitiesData = [0.3,-0.8];            
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

        
        function createOrientation(obj)
            m = obj.mesh;
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

        function plotOrientation(obj)
            figure()
            x = obj.mesh.coord(:,1);
            y = obj.mesh.coord(:,2);
            t  = obj.orientation;
            ct = cos(t(:,1));
            st = sin(t(:,1));
            quiver(x,y,ct,st)
        end        

         function s = createLevelSetCellParams(obj)        
            I        = ones(size(obj.mesh.coord,1),1);
            s.type   = 'rectangleInclusion';
            s.widthH = obj.widthH*I;
            s.widthV = obj.widthW*I;
            s.ndim   = 2;
         end       

        function dehomogenize(obj)
            s.nCells             = obj.nCells;
            s.cellLevelSetParams = obj.createLevelSetCellParams();
            s.mesh               = obj.mesh;
            s.theta              = atan2(obj.orientation(:,2),obj.orientation(:,1));
            d = Dehomogenizer(s);
            ls = d.compute();
            d.plot();
            obj.levelSet = ls;
        end         
        
    end
    
end