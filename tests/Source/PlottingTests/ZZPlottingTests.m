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

    methods (Test, TestTags = {'PlottingTests', 'Passed', 'Classic', 'Displacement'})

        function testDisplacement(testCase, coordTests)
            cd ../../../
            s.testName              = coordTests;
            s.meshType              = 'BOUNDARY';
            s.meshIncludeBoxContour = false;
            test = ZZTestPlotting(s);
            err = test.computeError();
            tol = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
            cd ./tests/Source/PlottingTests/
        end

    end
end