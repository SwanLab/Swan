classdef AcademicTests < matlab.unittest.TestCase

    properties (TestParameter)
        problem = {'AcademicTest0','AcademicTest1', 'AcademicTest2', 'AcademicTest3'...
            'AcademicTest4'}
    end

    methods (Test, TestTags = {'Academic'})

        function testFast(testCase, problem)
            run(problem);
            cParams.cost         = cost;
            cParams.constraint   = constraint;
            cParams.initialGuess = x0;
            cParams.settings     = s;
            cParams.printingPath = false;
            test                 = AcademicProblem(cParams);
            test.compute();
            obt = test.result;
            org = load(problem).result;
            err = abs(obt.fun.fValues-org);
            tol = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
        end

    end

end
