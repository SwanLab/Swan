classdef DehomogenizingSingularitiesTest < handle

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
        alpha
    end

    methods (Access = public)

        function obj = DehomogenizingSingularitiesTest(cParams)
            obj.init(cParams);
            %mSize = 0.02;
            mSize = linspace(0.02,0.03,2);
            %mSize = linspace(0.030526315789474,0.038421052631579,0.034137931034483,0.034482758620690,0.036551724137931);
         %   mSize = linspace(0.04,0.05,30);%linspace(0.042,0.04,2);%0.09;%0.0221;%0.09;%0.0221;%0.09;%0.0221;%0.0521 %0.0221;0.0921
            for iMesh = 1:length(mSize)
                iMesh
                obj.meshSize = mSize(iMesh);
                obj.createMesh();
                obj.createAngle();
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
            obj.meshSize = 0.06;%0.0221;%0.09;%0.0221;%0.09;%0.0221;%0.0521 %0.0221;0.0921
            obj.nCells   = 30;%linspace(10,62,40); %[60 62];
           % 
           %  obj.xmin = 0.5;
           % obj.xmax = 2.0;
           % obj.ymin = 0.25;
           % obj.ymax = 1.75;
           
             obj.xmin = -0.2;%.2;
             obj.xmax = 2.0;
             obj.ymin = 0.25;
             obj.ymax = 1.75;

             % obj.xmin = 0.235;
             % obj.xmax = 1.75;
             % obj.ymin = -0.6;
             % obj.ymax = 0.705; 

            obj.singularitiesData = [0.3,-0.8];%[0.32,-0.8];
            obj.widthH = 0.89;
            obj.widthW = 0.89;
        end

        function createMesh(obj)
%            [p,b,t,nv,nbe,nt,labels]=ffreadmesh('export_mesh.msh');

            % t = t';
            % p = p';
            % s.connec = t(:,1:3);
            % s.coord  = p(:,1:2);
            % m = Mesh(s);
            % obj.mesh = m;
            % 

            h = obj.meshSize;
            xv = obj.xmin:h:obj.xmax;
            yv = obj.ymin:h:obj.ymax;
            [X,Y] = meshgrid(xv,yv);
            s.coord(:,1) = X(:);
            s.coord(:,2) = Y(:);
              [F,V] = mesh2tri(X,Y,zeros(size(X)),'f');
              s.coord  = V(:,1:2);
              s.connec = F;

%            s.connec = delaunay(s.coord);
            m = Mesh.create(s);
            obj.mesh = m;
            
        %    s.filename = 'SingMeshLeft';
        %    s.type     = 'GiD';
        %    m.print(s);

      %     m.exportSTL('SingMesh')
      %     mesh = msh('SingMesh.stl');
      %    mshWriteMsh('/home/alex/Desktop/PerleSingularities/SingMeshLeft.msh',mesh)
            

       %     s.filename= 'SingMesh';
       %     s.type = 'GiD';
       %     m.print(s)            

        end

        function createAngle(obj)
            m = obj.mesh;
            s1 = obj.singularitiesData(:,1);
            s2 = obj.singularitiesData(:,2);
            x1 = m.coord(:,1);
            x2 = m.coord(:,2);
            v(:,1) = cos(pi*(x1 + s1*x2));
            v(:,2) = cos(pi*(x2 + s2*x1));         
            beta = atan2(v(:,2),v(:,1));
            alpha = beta/2;
            s.fValues = alpha;
            s.mesh    = obj.mesh;
            s.order   = 'P1';
            aF = LagrangianFunction(s);
            obj.alpha = aF;
        end


        function createOrientation(obj)
            al = obj.alpha.fValues;
            a = [cos(al), sin(al)];
            s.fValues = a;
            s.mesh    = obj.mesh;
            s.order   = 'P1';
            aF = LagrangianFunction(s);            
            obj.orientation{1} = aF;
            a = [-sin(al), cos(al)];
            s.fValues = a;
            s.mesh    = obj.mesh;
            aF = LagrangianFunction(s);            
            obj.orientation{2} = aF;
        end

        function plotOrientation(obj)
            figure()
            x = obj.mesh.coord(:,1);
            y = obj.mesh.coord(:,2);
            t  = obj.orientation.fValues;
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
            s.theta              = obj.orientation;%.project('P0');%atan2(obj.orientation(:,2),obj.orientation(:,1));
            d = Dehomogenizer(s);
            ls = d.compute();
            d.plot();
            obj.levelSet = ls;
        end

    end

end
