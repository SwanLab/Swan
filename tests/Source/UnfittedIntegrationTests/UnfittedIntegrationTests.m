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
            a.fileName = filename;
            f = FemDataContainer(a);
            mesh = f.mesh;

            s.testName              = circleTests;
            s.analyticalValue       = 2*pi;
            s.meshType              = 'BOUNDARY';
            s.meshIncludeBoxContour = false;
            test = UnfittedIntegrationTest(s);
            err = test.computeError();
            tol = 6e-2;
            testCase.verifyLessThanOrEqual(err, tol)
        end

    end

    methods (Test, TestTags = {'UnfittedIntegration', 'Classic', 'Area', 'Circle'})

        function testAreaCircle(testCase, circleTests)
            testCase.fixFolder();

            run(circleTests)
            a.fileName = filename;
            f = FemDataContainer(a);
            mesh = f.mesh;

            s.testName              = circleTests;
            s.analyticalValue       = pi;
            s.meshType              = 'INTERIOR';
            s.meshIncludeBoxContour = false;
            test = UnfittedIntegrationTest(s);
            err = test.computeError();
            tol = 6e-2;
            testCase.verifyLessThanOrEqual(err, tol)
        end

    end

    methods (Test, TestTags = {'UnfittedIntegration', 'Classic', 'Perimeter', 'Rectangle'})

        function testRectangle(testCase, rectangleTests)
            testCase.fixFolder();

            run(rectangleTests)
            a.fileName = filename;
            f = FemDataContainer(a);
            mesh = f.mesh;

            s.testName              = rectangleTests;
            s.analyticalValue       = 6;
            s.meshType              = 'BOUNDARY';
            s.meshIncludeBoxContour = true;
            test = UnfittedIntegrationTest(s);
            err = test.computeError();
            tol = 6e-2;
            testCase.verifyLessThanOrEqual(err, tol)
        end

    end

    methods (Test, TestTags = {'UnfittedIntegration', 'Classic', 'Surface', 'Sphere'})

        function testSurfaceSphere(testCase, sphereTests)
            testCase.fixFolder();

            run(sphereTests)
            a.fileName = filename;
            f = FemDataContainer(a);
            mesh = f.mesh;

            s.testName              = sphereTests;
            s.analyticalValue       = 4*pi;
            s.meshType              = 'BOUNDARY';
            s.meshIncludeBoxContour = false;
            test = UnfittedIntegrationTest(s);
            err = test.computeError();
            tol = 6e-2;
            testCase.verifyLessThanOrEqual(err, tol)
        end

    end

    methods (Test, TestTags = {'UnfittedIntegration', 'Classic', 'Volume', 'Sphere'})

        function testVolumeSphere(testCase, sphereTests)
            testCase.fixFolder();

            run(sphereTests)
            a.fileName = filename;
            f = FemDataContainer(a);
            mesh = f.mesh;

            s.testName              = sphereTests;
            s.analyticalValue       = (4/3)*pi;
            s.meshType              = 'INTERIOR';
            s.meshIncludeBoxContour = false;
            test = UnfittedIntegrationTest(s);
            err = test.computeError();
            tol = 6e-2;
            testCase.verifyLessThanOrEqual(err, tol)
        end

    end

    methods (Test, TestTags = {'UnfittedIntegration', 'Classic', 'Surface', 'Cylinder'})

        function testSurfaceCylinder(testCase, cylinderTests)
            testCase.fixFolder();
            run(cylinderTests)
            a.fileName = filename;
            f = FemDataContainer(a);
            mesh = f.mesh;
            gPar.type               = 'Cylinder';
            gPar.radius             = 1;
            gPar.xCoorCenter        = 1;
            gPar.yCoorCenter        = 1;
            g                       = GeometricalFunction(gPar);
            phiFun                  = g.computeLevelSetFunction(mesh);
            phi                     = phiFun.fValues;
            s.analyticalValue       = pi*2+2*pi*2;
            s.meshType              = 'BOUNDARY';
            s.mesh                  = mesh;
            s.levelSet              = phi;
            test = UnfittedIntegrationTest(s);
            err = test.computeError();
            tol = 6e-2;
            testCase.verifyLessThanOrEqual(err, tol)
        end

    end

    methods (Test, TestTags = {'UnfittedIntegration', 'Nou', 'Volume', 'Cylinder'})

        function testVolumeCylinder(testCase, cylinderTests)
            testCase.fixFolder();
            run(cylinderTests)
            a.fileName = filename;
            f = FemDataContainer(a);
            mesh = f.mesh;
            gPar.type               = 'Cylinder';
            gPar.radius             = 1;
            gPar.xCoorCenter        = 1;
            gPar.yCoorCenter        = 1;
            g                       = GeometricalFunction(gPar);
            phiFun                  = g.computeLevelSetFunction(mesh);
            phi                     = phiFun.fValues;
            s.analyticalValue       = pi*2;
            s.meshType              = 'INTERIOR';
            s.mesh                  = mesh;
            s.levelSet              = phi;
            test = UnfittedIntegrationTest(s);
            err = test.computeError();
            tol = 6e-2;
            testCase.verifyLessThanOrEqual(err, tol)
        end

    end

%     methods (Test, TestTags = {'UnfittedIntegration', 'Classic', 'Line', 'halfCylinder'})
% 
%         function testLinesHalfCylinder(testCase, halfCylinderTests)
%             testCase.fixFolder();
%             s.testName              = halfCylinderTests;
%             s.analyticalValue       = pi*2 + 2*pi*2;
%             s.meshType              = 'BOUNDARY';
%             s.meshIncludeBoxContour = true; % irrelevant
%             test = UnfittedIntegrationTest(s);
%             err = test.computeError();
%             tol = 6e-2;
%             testCase.verifyLessThanOrEqual(err, tol)
%         end
% 
%     end
% 
%     methods (Test, TestTags = {'UnfittedIntegration', 'Classic', 'Surface', 'halfCylinder'})
% 
%         function testSurfaceHalfCylinder(testCase, halfCylinderTests)
%             testCase.fixFolder();
%             s.testName              = halfCylinderTests;
%             s.analyticalValue       = pi*2 + 2*pi*2;
%             s.meshType              = 'BOUNDARY';
%             s.meshIncludeBoxContour = true; % irrelevant
%             test = UnfittedIntegrationTest(s);
%             err = test.computeError();
%             tol = 6e-2;
%             testCase.verifyLessThanOrEqual(err, tol)
%         end
% 
%     end
% 
%     methods (Test, TestTags = {'UnfittedIntegration', 'Nou', 'Volume', 'halfCylinder'})
% 
%         function testVolumeHalfCylinder(testCase, halfCylinderTests)
%             testCase.fixFolder();
%             s.testName              = halfCylinderTests;
%             s.analyticalValue       = pi*2;
%             s.meshType              = 'INTERIOR';
%             s.meshIncludeBoxContour = false; % irrelevant
%             test = UnfittedIntegrationTest(s);
%             err = test.computeError();
%             tol = 6e-2;
%             testCase.verifyLessThanOrEqual(err, tol)
%         end
% 
%     end

    methods (Access = private)

        function fixFolder(testCase)
            import matlab.unittest.fixtures.CurrentFolderFixture
            changeToFolder = '../../../';
            testCase.applyFixture(CurrentFolderFixture(changeToFolder));
        end
    end

end