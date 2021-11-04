classdef UnfittedIntegrationTests < handle & matlab.unittest.TestCase

    properties (TestParameter)
        circleTests = {'test_circle_triangle', 'test_circle_quadrilateral'}
        rectangleTests = {'test_rectangle_triangle', 'test_rectangle_quadrilateral'}
        sphereTests = {'test_sphere_hexahedra', 'test_sphere_tetrahedra'}
        cylinderTests = {'test_cylinder_tetrahedra', 'test_cylinder_hexahedra'}
    end

    methods (Test, TestTags = {'UnfittedIntegration', 'Classic', 'Perimeter', 'Circle'})
        function testPerimeterCircle(testCase, circleTests)
            % Changes in directories are used to fix the current directory
            % so that Unit Tests can pass. This should not be an issue when
            % OldSettings disappear.
            cd ../../../
            s.testName              = circleTests;
            s.analyticalValue       = 2*pi;
            s.meshType              = 'BOUNDARY';
            s.meshIncludeBoxContour = false;
            test = UnfittedIntegrationTest(s);
            err = test.computeError();
            tol = 6e-2;
            testCase.verifyLessThanOrEqual(err, tol)
            cd ./tests/Source/UnfittedIntegrationTests/
        end

    end

    methods (Test, TestTags = {'UnfittedIntegration', 'Classic', 'Area', 'Circle'})

        function testAreaCircle(testCase, circleTests)
            cd ../../../
            s.testName              = circleTests;
            s.analyticalValue       = pi;
            s.meshType              = 'INTERIOR';
            s.meshIncludeBoxContour = false;
            test = UnfittedIntegrationTest(s);
            err = test.computeError();
            tol = 6e-2;
            testCase.verifyLessThanOrEqual(err, tol)
            cd ./tests/Source/UnfittedIntegrationTests/
        end

    end

    methods (Test, TestTags = {'UnfittedIntegration', 'Classic', 'Perimeter', 'Rectanngle'})

        function testRectangle(testCase, rectangleTests)
            cd ../../../
            s.testName              = rectangleTests;
            s.analyticalValue       = 6;
            s.meshType              = 'BOUNDARY';
            s.meshIncludeBoxContour = true;
            test = UnfittedIntegrationTest(s);
            err = test.computeError();
            tol = 6e-2;
            testCase.verifyLessThanOrEqual(err, tol)
            cd ./tests/Source/UnfittedIntegrationTests/
        end

    end

    methods (Test, TestTags = {'UnfittedIntegration', 'Classic', 'Surface', 'Sphere'})

        function testSurfaceSphere(testCase, sphereTests)
            cd ../../../
            s.testName              = sphereTests;
            s.analyticalValue       = 4*pi;
            s.meshType              = 'BOUNDARY';
            s.meshIncludeBoxContour = false;
            test = UnfittedIntegrationTest(s);
            err = test.computeError();
            tol = 6e-2;
            testCase.verifyLessThanOrEqual(err, tol)
            cd ./tests/Source/UnfittedIntegrationTests/
        end

    end

    methods (Test, TestTags = {'UnfittedIntegration', 'Classic', 'Volume', 'Sphere'})

        function testVolumeSphere(testCase, sphereTests)
            cd ../../../
            s.testName              = sphereTests;
            s.analyticalValue       = (4/3)*pi;
            s.meshType              = 'INTERIOR';
            s.meshIncludeBoxContour = false;
            test = UnfittedIntegrationTest(s);
            err = test.computeError();
            tol = 6e-2;
            testCase.verifyLessThanOrEqual(err, tol)
            cd ./tests/Source/UnfittedIntegrationTests/
        end

    end

    methods (Test, TestTags = {'UnfittedIntegration', 'Classic', 'Surface', 'Cylinder'})

        function testSurfaceCylinder(testCase, cylinderTests)
            cd ../../../
            s.testName              = cylinderTests;
            s.analyticalValue       = pi*2 + 2*pi*2;
            s.meshType              = 'BOUNDARY';
            s.meshIncludeBoxContour = true; % irrelevant
            test = UnfittedIntegrationTest(s);
            err = test.computeError();
            tol = 6e-2;
            testCase.verifyLessThanOrEqual(err, tol)
            cd ./tests/Source/UnfittedIntegrationTests/
        end

    end

    methods (Test, TestTags = {'UnfittedIntegration', 'Nou', 'Volume', 'Cylinder'})

        function testVolumeCylinder(testCase, cylinderTests)
            cd ../../../
            s.testName              = cylinderTests;
            s.analyticalValue       = pi*2;
            s.meshType              = 'INTERIOR';
            s.meshIncludeBoxContour = false; % irrelevant
            test = UnfittedIntegrationTest(s);
            err = test.computeError();
            tol = 6e-2;
            testCase.verifyLessThanOrEqual(err, tol)
            cd ./tests/Source/UnfittedIntegrationTests/
        end

    end

end