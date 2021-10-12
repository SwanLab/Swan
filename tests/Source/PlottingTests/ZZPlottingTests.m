classdef ZZPlottingTests < handle & matlab.unittest.TestCase

    properties (TestParameter)
        coordTests = {'test_triangleToyUntittedExample', 'test_quadToyUntittedExample'} %no file
        plottingTests = {'test_circumference_quadrilateral', 'test_circumference_triangle', ...
            'test_circle_triangle', 'test_circle_quadrilateral', ...
            'test_sphere_tetrahedra', 'test_sphere_hexahedra'}
        compositeTests = {'test_rectangle_triangle_plot', 'test_rectangle_quadrilateral_plot', ...
            'test_smoothRectangle_triangle', 'test_smoothRectangle_quadrilateral', ...
            'test_cylinder_tetrahedra', 'test_cylinder_hexahedra'}
        boletTests = {'testPlotLargeCylinderTethaedra'}
    end

    methods (Test, TestTags = {'PlottingTests', 'Toy', 'Nou'})

        function testTriangleToy(testCase, coordTests)
            cd ../../../
            s.testName              = coordTests;
            s.meshType              = 'BOUNDARY';
            s.meshIncludeBoxContour = false;
            obj.coord  = [0 0;1 0;1 1;0 1;2 0;2 1;0 2;1 2;2 2; 0.5 0.5; 1.5 0.5; 0.5 1.5; 1.5 1.5];
            obj.connec = [1 2 10; 2 3 10; 10 3 4; 10 4 1; 2 11 3; 2 5 11; 5 6 11; 11 6 3; 3 8 12; 4 3 12; 12 8 7; 12 7 4; 3 6 13; 6 9 13; 13 9 8; 3 13 8];
            obj.levelSet = [-0.05 0.2 -0.5 0.1 0.1 -1 1 -0.2 -0.5 -0.05 -0.05 0.05 -0.5]';
            test = ZZTestPlotting(s);
            err = test.computeError();
            tol = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
            cd ./tests/Source/PlottingTests/
        end

    end

%     methods (Test, TestTags = {'PlottingTests', 'Passed', 'Classic', 'Displacement'})
% 
%         function testDisplacement(testCase, coordTests)
%             cd ../../../
%             s.testName              = coordTests;
%             s.meshType              = 'BOUNDARY';
%             s.meshIncludeBoxContour = false;
%             test = ZZTestPlotting(s);
%             err = test.computeError();
%             tol = 1e-6;
%             testCase.verifyLessThanOrEqual(err, tol)
%             cd ./tests/Source/PlottingTests/
%         end
% 
%     end
end