classdef PlottingTests < handle & matlab.unittest.TestCase

    properties (TestParameter)
        % coordTests = {'test_triangleToyUntittedExample', 'test_quadToyUntittedExample'} %no file
        plottingTests = {'test_circumference_quadrilateral', 'test_circumference_triangle', ...
            'test_circle_triangle', 'test_circle_quadrilateral', ...
            'test_sphere_tetrahedra', 'test_sphere_hexahedra'}
        compositeTests = {'test_rectangle_triangle', 'test_rectangle_quadrilateral', ...
            'test_smoothRectangle_triangle', 'test_smoothRectangle_quadrilateral', ...
            'test_cylinder_tetrahedra', 'test_cylinder_hexahedra'}
        % cylinderTests = {'testPlotLargeCylinderTethaedra'}
    end


%     methods (Test, TestTags = {'PlottingTests', 'ToPass'})
% 
%         function testsThatFail(testCase, plottingTests)
%             testCase.fixFolder();
%             s.testName              = plottingTests;
%             s.meshType              = 'BOUNDARY';
%             s.meshIncludeBoxContour = false;
%             test = TestPlotting(s);
%             passed = test.computePassed();
%             verifyTrue(testCase, passed)
%         end
% 
%     end

    methods (Test, TestTags = {'PlottingTests', 'Toy'})

        function testTriangleToy(testCase)
            testCase.fixFolder();
            s.testName = 'test_triangleToyUntittedExample';
            s.coord  = [0 0;1 0;1 1;0 1;2 0;2 1;0 2;1 2;2 2; 0.5 0.5;1.5 0.5; 0.5 1.5; 1.5 1.5];
            s.connec = [1 2 10; 2 3 10; 10 3 4; 10 4 1; 2 11 3; 2 5 11; 5 6 11; 11 6 3; 3 8 12;
                        4 3 12; 12 8 7; 12 7 4; 3 6 13; 6 9 13; 13 9 8; 3 13 8];
            s.levelSet = [-0.05 0.2 -0.5 0.1 0.1 -1 1 -0.2 -0.5 -0.05 -0.05 0.05 -0.5]';
            test = PlottingToyUnfittedExample(s);
            passed = test.hasPassed();
            verifyTrue(testCase, passed)
        end

        function testQuadToy(testCase)
            testCase.fixFolder();
            s.testName = 'test_quadToyUntittedExample';
            s.coord  = [0 0;1 0;1 1;0 1;2 0;2 1;0 2;1 2;2 2];
            s.connec = [1 2 3 4; 2 5 6 3; 4 3 8 7; 3 6 9 8];
            s.levelSet = [-0.05 0.2 -0.5 0.1 0.1 -1 1 -0.2 -0.5]';
            test = PlottingToyUnfittedExample(s);
            passed = test.hasPassed();
            verifyTrue(testCase, passed)
        end

    end

    methods (Test, TestTags = {'PlottingTests', 'FileBased'})

        function testPlottingNormal(testCase, plottingTests)
            testCase.fixFolder();
            s.testName              = plottingTests;
            s.meshType              = 'BOUNDARY';
            s.meshIncludeBoxContour = false;
            test = TestPlotting(s);
            passed = test.computePassed();
            verifyTrue(testCase, passed)
        end

        function testComposite(testCase, compositeTests)
            testCase.fixFolder();
            s.testName              = compositeTests;
            s.meshType              = 'BOUNDARY';
            s.meshIncludeBoxContour = true;
            test = TestPlotting(s);
            passed = true;
            verifyTrue(testCase, passed)
        end

    end

    % methods (Test, TestTags = {'PlottingTests', 'Cylinder'})
    % 
    %     function testCylinder(testCase)
    %         test = TestPlotLargeCylinderTethaedra();
    %         passed = true;
    %         verifyTrue(testCase, passed)
    %     end
    % 
    % end

    methods (Access = private)
        
        function fixFolder(testCase)
            import matlab.unittest.fixtures.CurrentFolderFixture
            changeToFolder = '../../../';
            testCase.applyFixture(CurrentFolderFixture(changeToFolder));
        end
    end

end