classdef UnfittedIntegrationTests < handle & matlab.unittest.TestCase

    properties (TestParameter)
        circleTests = {'test_circle_triangle', 'test_circle_quadrilateral'}
        rectangleTests = {'test_rectangle_triangle', 'test_rectangle_quadrilateral'}
        sphereTests = {'test_sphere_hexahedra', 'test_sphere_tetrahedra'}
        cylinderTests = {'test_cylinder_tetrahedra', 'test_cylinder_hexahedra'}
    end

    methods (Test, TestTags = {'UnfittedIntegration', 'Classic', 'Perimeter', 'Circle'})
        function testPerimeterCircle(testCase, circleTests)
            testCase.fixFolder();
            run(circleTests)
            gPar.type         = 'Circle';
            gPar.radius       = 1;
            gPar.xCoorCenter  = 1;
            gPar.yCoorCenter  = 1;
            g                 = GeometricalFunction(gPar);
            phiFun            = g.computeLevelSetFunction(mMesh);
            phi               = phiFun.fValues;
            s.analyticalValue = 2*pi;
            s.meshType        = 'BOUNDARY';
            s.mesh            = mMesh;
            s.levelSet        = phi;
            test = UnfittedIntegrationTest(s);
            err  = test.computeError();
            tol  = 6e-2;
            testCase.verifyLessThanOrEqual(err, tol)
        end

    end

    methods (Test, TestTags = {'UnfittedIntegration', 'Classic', 'Area', 'Circle'})

        function testAreaCircle(testCase, circleTests)
            testCase.fixFolder();
            run(circleTests)
            gPar.type         = 'Circle';
            gPar.radius       = 1;
            gPar.xCoorCenter  = 1;
            gPar.yCoorCenter  = 1;
            g                 = GeometricalFunction(gPar);
            phiFun            = g.computeLevelSetFunction(mMesh);
            phi               = phiFun.fValues;
            s.analyticalValue = pi;
            s.meshType        = 'INTERIOR';
            s.mesh            = mMesh;
            s.levelSet        = phi;
            test = UnfittedIntegrationTest(s);
            err  = test.computeError();
            tol  = 6e-2;
            testCase.verifyLessThanOrEqual(err, tol)
        end

    end

    methods (Test, TestTags = {'UnfittedIntegration', 'Classic', 'Perimeter', 'Rectangle'})

        function testRectangle(testCase, rectangleTests)
            testCase.fixFolder();
            run(rectangleTests)
            gPar.type         = 'Rectangle';
            gPar.xSide        = 1;
            gPar.ySide        = 2;
            gPar.xCoorCenter  = 1;
            gPar.yCoorCenter  = 1;
            g                 = GeometricalFunction(gPar);
            phiFun            = g.computeLevelSetFunction(mMesh);
            phi               = phiFun.fValues;
            s.analyticalValue = 6;
            s.meshType        = 'BOUNDARY';
            s.mesh            = mMesh;
            s.levelSet        = phi;
            test = UnfittedIntegrationTest(s);
            err  = test.computeError();
            tol  = 6e-2;
            testCase.verifyLessThanOrEqual(err, tol)
        end

    end

    methods (Test, TestTags = {'UnfittedIntegration', 'Classic', 'Surface', 'Sphere'})

        function testSurfaceSphere(testCase, sphereTests)
            testCase.fixFolder();
            run(sphereTests)
            gPar.type         = 'Sphere';
            gPar.radius       = 1;
            gPar.xCoorCenter  = 1;
            gPar.yCoorCenter  = 1;
            gPar.zCoorCenter  = 1;
            g                 = GeometricalFunction(gPar);
            phiFun            = g.computeLevelSetFunction(mMesh);
            phi               = phiFun.fValues;
            s.analyticalValue = 4*pi;
            s.meshType        = 'BOUNDARY';
            s.mesh            = mMesh;
            s.levelSet        = phi;
            test = UnfittedIntegrationTest(s);
            err  = test.computeError();
            tol  = 6e-2;
            testCase.verifyLessThanOrEqual(err, tol)
        end

    end

    methods (Test, TestTags = {'UnfittedIntegration', 'Classic', 'Volume', 'Sphere'})

        function testVolumeSphere(testCase, sphereTests)
            testCase.fixFolder();
            run(sphereTests)
            gPar.type         = 'Sphere';
            gPar.radius       = 1;
            gPar.xCoorCenter  = 1;
            gPar.yCoorCenter  = 1;
            gPar.zCoorCenter  = 1;
            g                 = GeometricalFunction(gPar);
            phiFun            = g.computeLevelSetFunction(mMesh);
            phi               = phiFun.fValues;
            s.analyticalValue = (4/3)*pi;
            s.meshType        = 'INTERIOR';
            s.mesh            = mMesh;
            s.levelSet        = phi;
            test = UnfittedIntegrationTest(s);
            err  = test.computeError();
            tol  = 6e-2;
            testCase.verifyLessThanOrEqual(err, tol)
        end

    end

    methods (Test, TestTags = {'UnfittedIntegration', 'Classic', 'Surface', 'Cylinder'})

        function testSurfaceCylinder(testCase, cylinderTests)
            testCase.fixFolder();
            run(cylinderTests)
            gPar.type         = 'Cylinder';
            gPar.radius       = 1;
            gPar.xCoorCenter  = 1;
            gPar.yCoorCenter  = 1;
            g                 = GeometricalFunction(gPar);
            phiFun            = g.computeLevelSetFunction(mMesh);
            phi               = phiFun.fValues;
            s.analyticalValue = pi*2+2*pi*2;
            s.meshType        = 'BOUNDARY';
            s.mesh            = mMesh;
            s.levelSet        = phi;
            test = UnfittedIntegrationTest(s);
            err  = test.computeError();
            tol  = 6e-2;
            testCase.verifyLessThanOrEqual(err, tol)
        end

    end

    methods (Test, TestTags = {'UnfittedIntegration', 'Nou', 'Volume', 'Cylinder'})

        function testVolumeCylinder(testCase, cylinderTests)
            testCase.fixFolder();
            run(cylinderTests)
            gPar.type         = 'Cylinder';
            gPar.radius       = 1;
            gPar.xCoorCenter  = 1;
            gPar.yCoorCenter  = 1;
            g                 = GeometricalFunction(gPar);
            phiFun            = g.computeLevelSetFunction(mMesh);
            phi               = phiFun.fValues;
            s.analyticalValue = pi*2;
            s.meshType        = 'INTERIOR';
            s.mesh            = mMesh;
            s.levelSet        = phi;
            test = UnfittedIntegrationTest(s);
            err  = test.computeError();
            tol  = 6e-2;
            testCase.verifyLessThanOrEqual(err, tol)
        end

    end

    methods (Test, TestTags = {'UnfittedIntegration', 'Classic', 'Line', 'halfCylinder'})

        function testLinesHalfCylinder(testCase, cylinderTests)
            testCase.fixFolder();
            run(cylinderTests)
            gPar.type         = 'Cylinder';
            gPar.radius       = 0.8;
            gPar.xCoorCenter  = 2;
            gPar.yCoorCenter  = 1;
            g                 = GeometricalFunction(gPar);
            phiFun            = g.computeLevelSetFunction(mMesh);
            phi               = phiFun.fValues;
            s.analyticalValue = 2*pi*0.8+2+2+1.6+1.6;
            s.meshType        = 'SUBBOUNDARY';
            s.mesh            = mMesh;
            s.levelSet        = phi;
            test = UnfittedIntegrationTest(s);
            err  = test.computeError();
            tol  = 6e-2;
            testCase.verifyLessThanOrEqual(err, tol)
        end

    end

    methods (Test, TestTags = {'UnfittedIntegration', 'Classic', 'Surface', 'halfCylinder'})

        function testSurfaceHalfCylinder(testCase, cylinderTests)
            testCase.fixFolder();
            run(cylinderTests)
            gPar.type         = 'Cylinder';
            gPar.radius       = 0.8;
            gPar.xCoorCenter  = 2;
            gPar.yCoorCenter  = 1;
            g                 = GeometricalFunction(gPar);
            phiFun            = g.computeLevelSetFunction(mMesh);
            phi               = phiFun.fValues;
            s.analyticalValue = pi*0.8^2+2*1.6+2*pi*0.8;
            s.meshType        = 'BOUNDARY';
            s.mesh            = mMesh;
            s.levelSet        = phi;
            test = UnfittedIntegrationTest(s);
            err  = test.computeError();
            tol  = 6e-2;
            testCase.verifyLessThanOrEqual(err, tol)
        end

    end

    methods (Test, TestTags = {'UnfittedIntegration', 'Nou', 'Volume', 'halfCylinder'})

        function testVolumeHalfCylinder(testCase, cylinderTests)
            testCase.fixFolder();
            run(cylinderTests)
            gPar.type         = 'Cylinder';
            gPar.radius       = 0.8;
            gPar.xCoorCenter  = 2;
            gPar.yCoorCenter  = 1;
            g                 = GeometricalFunction(gPar);
            phiFun            = g.computeLevelSetFunction(mMesh);
            phi               = phiFun.fValues;
            s.analyticalValue = pi*0.8^2;
            s.meshType        = 'INTERIOR';
            s.mesh            = mMesh;
            s.levelSet        = phi;
            test = UnfittedIntegrationTest(s);
            err  = test.computeError();
            tol  = 6e-2;
            testCase.verifyLessThanOrEqual(err, tol)
        end

    end

    methods (Access = private)

        function fixFolder(testCase)
            import matlab.unittest.fixtures.CurrentFolderFixture
            changeToFolder = '../../../';
            testCase.applyFixture(CurrentFolderFixture(changeToFolder));
        end
    end

end