classdef AcademicTests < matlab.unittest.TestCase

    properties (TestParameter)
        problem = {'AcademicTest0','AcademicTest1', 'AcademicTest2', 'AcademicTest3'...
            'AcademicTest4'}
    end

    methods (Test, TestTags = {'Academic'})

        function testFast(testCase, problem)
            run(problem);
            cParams.cost           = cost;
            cParams.constraint     = constraint;
            cParams.constraint.nSF = nConstr;
            cParams.initialGuess   = x0;
            cParams.settings       = s;
            test                   = AcademicProblem(cParams);
            test.compute();
            obt = test.result;
            org = load(problem).result;
            err = abs(obt-org);
            tol = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
        end

    end

end
