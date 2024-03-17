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
            mSize = linspace(0.08,0.022,2);%0.09;%0.0221;%0.09;%0.0221;%0.09;%0.0221;%0.0521 %0.0221;0.0921
            for iMesh = 1:length(mSize)
                obj.meshSize = mSize(iMesh);
                obj.createMesh();
                obj.createOrientation();
                obj.dehomogenize();
            end
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
            obj.meshSize = 0.065;%0.09;%0.0221;%0.09;%0.0221;%0.09;%0.0221;%0.0521 %0.0221;0.0921
            obj.nCells   = 30;%[60 62];%linspace(60,62,40);%45;   %45        
            obj.xmin = 0.7;
            obj.xmax = 2.3;
            obj.ymin = 0.1;
            obj.ymax = 1.7;
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
           m = Mesh.create(s);
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
            al = beta/2;
            a = [cos(al), sin(al)];
            s.fValues = a;
            s.mesh    = obj.mesh;
            s.order   = 'P1';
            aF = LagrangianFunction(s);            
            obj.orientation{1} = aF;%.project('P0');
            a = [-sin(al), cos(al)];
            s.fValues = a;
            s.mesh    = obj.mesh;
            aF = LagrangianFunction(s);            
            obj.orientation{2} = aF;%.project('P0');
           % obj.orientation = a;
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
            s.theta              = obj.orientation;
            d = Dehomogenizer(s);
            ls = d.compute();
            d.plot();
            obj.levelSet = ls;
        end         
        
    end
    
end