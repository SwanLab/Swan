classdef TopOptTestOld < handle & matlab.unittest.TestCase

    properties (TestParameter)
        compTestsDehomog = {'test_M1M2', 'test_StressM1M2'}
        compTestsToPass = {'test_bridge', 'test_bridge2', ...
            'test_cantilever', 'test_cantilever2', 'test_cantilever3', ...
            'test_micro', 'test_micro2', ...
            'testDualNestedInPrimal_WithProjectedGradient', ...
            'testDualNestedInPrimal_WithSlerp', ...
            'test_interiorPerimeter'
            }
        fastDisp = {
            'test_cantilever', 'test_cantilever2', 'test_cantilever3', ...
            'test_micro', 'test_micro2', ...
            'testDualNestedInPrimal_WithProjectedGradient', ...
            'testDualNestedInPrimal_WithSlerp', ...
            'test_interiorPerimeter',...
            'test_anisotropy','test_anisotropy_interior','test_nullspace','test_gripping',...
            'test_interiorPerimeterPDErho','test_filterLump','test_cantilever_IPM'
            }
        macro = {'test_bridge2', ...
            'test_cantilever', 'test_cantilever2', 'test_cantilever3', ...
            'testDualNestedInPrimal_WithProjectedGradient', ...
            'testDualNestedInPrimal_WithSlerp', ...
            'test_interiorPerimeter'
            }
       micro = {'test_micro', 'test_micro2'}
%          micro = {'test_micro2'}
%         compTestsToPass = {'test_bridge'}
%         compTestsToPass = {'test_interiorPerimeter'}
        cantileverTests = {'test_cantilever2', 'test_cantilever3'}
        dimensions = {'2D', '3D'}
        vigdergauzTests = {'test_VigergauzMicroStructure', 'test_VigergauzMicroStructureFromStrain'}
        vigdergauzVolumes = {0.6, 0.75}
    end

    methods (Test, TestTags = {'Macro'})

        function testMacro(testCase, macro)
%             testCase.fixFolder();
            s.computerType    = 'TOPOPT';
            s.testName         = macro;
            s.variablesToStore = {'x'};
            test = PrecomputedVariableTest(s);
            err = test.computeError();
            tol = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
        end

    end

    methods (Test, TestTags = {'Micro'})

        function testMicro(testCase, micro)
%             testCase.fixFolder();
            s.computerType    = 'TOPOPT';
            s.testName         = micro;
            s.variablesToStore = {'x'};
            test = PrecomputedVariableTest(s);
            err = test.computeError();
            tol = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
        end

    end
    methods (Test, TestTags = {'TopOpt', 'Various', 'FastDisp'})

        function testFastDisplacement(testCase, fastDisp)
%             testCase.fixFolder();
            s.computerType    = 'TOPOPT';
            s.testName         = fastDisp;
            s.variablesToStore = {'x'};
            test = PrecomputedVariableTest(s);
            err = test.computeError();
            tol = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
        end

    end

    methods (Test, TestTags = {'Cantilever'})

        function testCantilever(testCase, cantileverTests)
%             testCase.fixFolder();
            s.computerType    = 'TOPOPT';
            s.testName         = cantileverTests;
            s.variablesToStore = {'x'};
            test = PrecomputedVariableTest(s);
            err = test.computeError();
            tol = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
        end

    end

    methods (Test, TestTags = {'TopOpt', 'Various', 'ToPass','Slow'})

        function testDisplacement(testCase, compTestsToPass)
%             testCase.fixFolder();
            s.computerType    = 'TOPOPT';
            s.testName         = compTestsToPass;
            s.variablesToStore = {'x'};
            test = PrecomputedVariableTest(s);
            err = test.computeError();
            tol = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
        end

    end

    methods (Test, TestTags = {'TopOpt', 'AnalyticVsRegularizedPerimeter', 'Displacement', 'Fast', 'Fixture'})

        function testAnalyticVsRegularizedPerimeter(testCase)
            test = TestAnalyticVsRegularizedPerimeter();
            err = test.computeError();
            tol = 5e-2;
            testCase.verifyLessThanOrEqual(err, tol)
        end

    end

    methods (Test, TestTags = {'TopOpt', 'SuperEllipseExponent', 'Fast'})

        function testSuperEllipseExponent(testCase)
            testCase.fixFolder();
            test = TestSuperEllipseExponent();
            err = test.computeError();
            tol = 1e-12;
            testCase.verifyLessThanOrEqual(err, tol)
        end

    end

    methods (Test, TestTags = {'TopOpt', 'ExplicitVsImplicit', 'Fast'})

        function testSimp(testCase, dimensions)
            testCase.fixFolder();
            test = SimplAllTestExplicitVsImplicit(dimensions);
            err = test.computeError();
            tol = 1e-12;
            testCase.verifyLessThanOrEqual(err, tol)
        end

    end

    methods (Test, TestTags = {'TopOpt', 'Vigdergauz', 'Micro', 'Medium'})

        function testVigdergauz(testCase, vigdergauzTests, vigdergauzVolumes)
%             testCase.fixFolder();
            s.testName = vigdergauzTests;
            s.volume = vigdergauzVolumes;
            test = TestVigdergauzMicroStructure(s);
            err = test.computeError();
            tol = 2*1e-1;
            testCase.verifyLessThanOrEqual(err, tol)
        end

    end

    methods (Access = private)
        
        function fixFolder(testCase)
            import matlab.unittest.fixtures.CurrentFolderFixture
            changeToFolder = '../../';
            testCase.applyFixture(CurrentFolderFixture(changeToFolder));
        end
    end

end
