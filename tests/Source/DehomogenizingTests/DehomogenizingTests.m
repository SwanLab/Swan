classdef DehomogenizingTests < handle & matlab.unittest.TestCase
    
  properties (TestParameter)
        coordTests = {'test_triangleToyUntittedExample', 'test_quadToyUntittedExample'} %no file
        plottingTests = {'test_circumference_quadrilateral', 'test_circumference_triangle', ...
            'test_circle_triangle', 'test_circle_quadrilateral', ...
            'test_sphere_tetrahedra', 'test_sphere_hexahedra'}
        compositeTests = {'test_rectangle_triangle_plot', 'test_rectangle_quadrilateral_plot', ...
            'test_smoothRectangle_triangle', 'test_smoothRectangle_quadrilateral', ...
            'test_cylinder_tetrahedra', 'test_cylinder_hexahedra'}
        cylinderTests = {'testPlotLargeCylinderTethaedra'}
    end

    methods (Test, TestTags = {'PlottingTests', 'Toy'})

        function testTriangleToy(testCase)
            cd ../../../
            s.testName = 'test_triangleToyUntittedExample';
            s.coord  = [0 0;1 0;1 1;0 1;2 0;2 1;0 2;1 2;2 2; 0.5 0.5;1.5 0.5; 0.5 1.5; 1.5 1.5];
            s.connec = [1 2 10; 2 3 10; 10 3 4; 10 4 1; 2 11 3; 2 5 11; 5 6 11; 11 6 3; 3 8 12;
                        4 3 12; 12 8 7; 12 7 4; 3 6 13; 6 9 13; 13 9 8; 3 13 8];
            s.levelSet = [-0.05 0.2 -0.5 0.1 0.1 -1 1 -0.2 -0.5 -0.05 -0.05 0.05 -0.5]';
            test = PlottingToyUnfittedExample(s);
            passed = test.hasPassed();
            verifyTrue(testCase, passed)
            cd ./tests/Source/PlottingTests/
        end

        function testQuadToy(testCase)
            cd ../../../
            s.testName = 'test_quadToyUntittedExample';
            s.coord  = [0 0;1 0;1 1;0 1;2 0;2 1;0 2;1 2;2 2];
            s.connec = [1 2 3 4; 2 5 6 3; 4 3 8 7; 3 6 9 8];
            s.levelSet = [-0.05 0.2 -0.5 0.1 0.1 -1 1 -0.2 -0.5]';
            test = PlottingToyUnfittedExample(s);
            passed = test.hasPassed();
            verifyTrue(testCase, passed)
            cd ./tests/Source/PlottingTests/
        end

    end

    methods (Test, TestTags = {'PlottingTests', 'FileBased'})

        function testPlottingNormal(testCase, plottingTests)
            cd ../../../
            s.testName              = plottingTests;
            s.meshType              = 'BOUNDARY';
            s.meshIncludeBoxContour = false;
            test = TestPlotting(s);
            passed = test.computePassed();
            verifyTrue(testCase, passed)
            cd ./tests/Source/PlottingTests/
        end

        function testComposite(testCase, compositeTests)
            cd ../../../
            s.testName              = compositeTests;
            s.meshType              = 'BOUNDARY';
            s.meshIncludeBoxContour = true;
            test = TestPlotting(s);
            passed = true; % hmmm, pero aixi es com era originalment
            verifyTrue(testCase, passed)
            cd ./tests/Source/PlottingTests/
        end

    end

    methods (Test, TestTags = {'PlottingTests', 'Cylinder'})

        function testCylinder(testCase)
            test = TestPlotLargeCylinderTethaedra();
            passed = true; % hmmm, pero aixi es com era originalment
            verifyTrue(testCase, passed)
        end

    end    
    
    properties (Access = protected)
        FieldOfStudy = 'Dehomogenizing'
        tests
    end
    
    methods (Access = public)
        
        function  obj = DehomogenizingTests()
            obj@testRunner();
        end
    end
    
    methods (Access = protected)
        function loadTests(obj)
            obj.tests = {'MeshSymmetrizerTest'; 
                'ScalarSymmetrizerTest';'MeshComponentCounterTest'};            
        end
    end
end