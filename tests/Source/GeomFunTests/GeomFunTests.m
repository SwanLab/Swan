classdef GeomFunTests < handle & matlab.unittest.TestCase

    properties (TestParameter)
        initialCases2D = {'Square','SmoothSquare','SquareInclusion','Rectangle',...
            'RectangleInclusion','SmoothRectangle','RectangleRotated',...
            'Circle','CircleInclusion','VerticalFiber','VerticalNFibers',...
            'HorizontalFiber','HorizontalInclusion','HorizontalNFibers',...
            'Holes','Full','Given','Vigdergauz','PeriodicAndOriented'}

        initialCases3D = {'Cylinder','Sphere','SphereInclusion','Holes',...
            'Full','Given'}
    end

    methods (Test, TestTags = {'GeomFun', 'Classic', 'LevelSet', 'Initial'})
        function testLevelSet2D(testCase, initialCases2D)
            filename = ['testGeomFun',initialCases2D,'2D'];
            m        = testCase.obtain2DTestMesh();
            s        = testCase.obtain2DTestsParameters(testCase,m,initialCases2D);
            s.type   = initialCases2D;
            g        = GeometricalFunction(s);
            lsFun    = g.computeLevelSetFunction(m);
            xNew     = lsFun.fValues;
            load(filename,'xRef');
            err      = norm(xNew-xRef)/norm(xRef);
            tol      = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
        end

        function testLevelSet3D(testCase, initialCases3D)
            filename = ['testGeomFun',initialCases3D,'3D'];
            m        = testCase.obtain3DTestMesh();
            s        = testCase.obtain3DTestsParameters();
            s.type   = initialCases3D;
            g        = GeometricalFunction(s);
            lsFun    = g.computeLevelSetFunction(m);
            xNew     = lsFun.fValues;
            load(filename,'xRef');
            err      = norm(xNew-xRef)/norm(xRef);
            tol      = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
        end
    end

    methods (Static, Access = private)
        function m = obtain2DTestMesh()
            x1       = linspace(-1,1,20);
            x2       = linspace(0,1,20);
            [xv,yv]  = meshgrid(x1,x2);
            [F,V]    = mesh2tri(xv,yv,zeros(size(xv)),'x');
            s.coord  = V(:,1:2);
            s.connec = F;
            m        = Mesh.create(s);
        end

        function m = obtain3DTestMesh()
            file       = 'Sphere_Hexahedra_Linear_Structured_8';
            a.fileName = file;
            s          = FemDataContainer(a);
            m          = s.mesh;
        end

        function s = obtain2DTestsParameters(obj,mesh,testType)
            switch testType
                case 'Vigdergauz'
                    s.vigdergauzSettings.type              = 'VolumeAndRatio';
                    s.vigdergauzSettings.volumeMicro       = 0.75;
                    s.vigdergauzSettings.superEllipseRatio = 2.5;
                case 'PeriodicAndOriented'
                    s = obj.obtainPeriodicOrientedSettings(mesh);
                otherwise
                    s.xCoorCenter  = 0;
                    s.yCoorCenter  = 0.5;
                    s.length       = 0.5;
                    s.pnorm        = 4;
                    s.xSide        = 1;
                    s.ySide        = 0.5;
                    s.omegaDeg     = 10;
                    s.radius       = 0.25;
                    s.width        = 0.15;
                    s.nFibers      = 3;
                    s.minxCoor     = -1;
                    s.maxxCoor     = 1;
                    s.minyCoor     = 0;
                    s.maxyCoor     = 1;
                    s.dim          = 2;
                    s.nHoles       = [4,3];
                    s.totalLengths = [2,1];
                    s.phases       = [0,0];
                    s.phiZero      = 0.4;
                    s.fHandle      = @(x) 0.25+atan(3*x(1,:,:).*x(2,:,:));
            end
        end

        function s = obtain3DTestsParameters()
            s.xCoorCenter  = 1;
            s.yCoorCenter  = 1;
            s.zCoorCenter  = 1;
            s.radius       = 0.5;
            s.dim          = 3;
            s.nHoles       = [4,3,3];
            s.totalLengths = [2,1,1];
            s.phases       = [0,0,0];
            s.phiZero      = 0.4;
            s.fHandle      = @(x) -1+atan(3*x(1,:,:).*x(2,:,:).*x(3,:,:));
        end

        function s = obtainPeriodicOrientedSettings(mesh)
            s1     = 0.3;
            s2     = -0.8;
            x1     = mesh.coord(:,1);
            x2     = mesh.coord(:,2);
            v(:,1) = cos(pi*(x1 + s1*x2));
            v(:,2) = cos(pi*(x2 + s2*x1));
            beta   = atan2(v(:,2),v(:,1));
            al     = beta/2;

            ss.fValues     = [cos(al), sin(al)];
            ss.mesh        = mesh;
            ss.order       = 'P1';
            aF             = LagrangianFunction(ss);
            orientation{1} = aF;
            ss.fValues     = [-sin(al), cos(al)];
            aF             = LagrangianFunction(ss);
            orientation{2} = aF;

            nCells = 30;
            L = mesh.computeCharacteristicLength();
            e = L/nCells;
            epsilons = (1 + 0.01*linspace(0,100,1))*e;

            sM.orientationP1 = orientation;
            sM.mesh          = mesh;
            
            I                          = ones(size(mesh.coord,1),1);
            s.orientationVectors       = OrientedMappingComputer(sM);
            s.epsilons                 = epsilons;
            s.mesh                     = mesh;
            s.cellLevelSetParams.type  = 'RectangleInclusion';
            s.cellLevelSetParams.xSide = 1*I;
            s.cellLevelSetParams.ySide = 0.5*I;
            s.cellLevelSetParams.xCoorCenter = 0;
            s.cellLevelSetParams.yCoorCenter = 0.5;
        end
    end

end