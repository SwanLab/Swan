classdef AcademicTests < matlab.unittest.TestCase

    properties (TestParameter)
        tests = {'AcademicTest1', 'AcademicTest2', 'AcademicTest3'}
    end

    methods (Test, TestTags = {'Academic'})

        function testTriangle(testCase, tests)
            s.filename = tests;
            test = AcademicProblem(s);
            obt = test.result;
            org = load(tests).result;
            err = abs(obt-org);
            tol = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
        end

    end

end
