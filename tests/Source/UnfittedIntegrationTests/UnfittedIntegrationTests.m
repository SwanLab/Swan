classdef UnfittedIntegrationTests < handle & matlab.unittest.TestCase

    properties (TestParameter)
        circleTests = {'test_circle_triangle', 'test_circle_quadrilateral'}
        rectangleTests = {'test_rectangle_triangle', 'test_rectangle_quadrilateral'}
        sphereTests = {'test_sphere_hexahedra', 'test_sphere_tetrahedra'}
        cylinderTests = {'test_cylinder_tetrahedra', 'test_cylinder_hexahedra'}
    end

    methods (Test, TestTags = {'UnfittedIntegration', 'Classic', 'Perimeter', 'Circle'})
        function testPerimeterCircle(testCase, circleTests)
            cd ../../../ % NOOOO, al segon test torna a saltar enrere
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
            cd ../../../ % NOOOO, al segon test torna a saltar enrere
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
            cd ../../../ % NOOOO, al segon test torna a saltar enrere
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
            cd ../../../ % NOOOO, al segon test torna a saltar enrere
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
            cd ../../../ % NOOOO, al segon test torna a saltar enrere
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
            cd ../../../ % NOOOO, al segon test torna a saltar enrere
            s.testName              = cylinderTests;
            s.analyticalValue       = pi*2 + 2*pi*2;
            s.meshType              = 'BOUNDARY';
            s.meshIncludeBoxContour = true; % no hi ha diferencia amb false
            test = UnfittedIntegrationTest(s);
            err = test.computeError();
            tol = 6e-2;
            testCase.verifyLessThanOrEqual(err, tol)
            cd ./tests/Source/UnfittedIntegrationTests/
        end

    end

    methods (Test, TestTags = {'UnfittedIntegration', 'Nou', 'Volume', 'Cylinder'})

        function testVolumeCylinder(testCase, cylinderTests)
            cd ../../../ % NOOOO, al segon test torna a saltar enrere
            s.testName              = cylinderTests;
            s.analyticalValue       = pi*2;
            s.meshType              = 'INTERIOR';
            s.meshIncludeBoxContour = false; % no hi ha diferencia amb true
            test = UnfittedIntegrationTest(s);
            err = test.computeError();
            tol = 6e-2;
            testCase.verifyLessThanOrEqual(err, tol)
            cd ./tests/Source/UnfittedIntegrationTests/
        end

    end

end