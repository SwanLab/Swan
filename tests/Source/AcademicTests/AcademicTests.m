classdef AcademicTests < matlab.unittest.TestCase

    properties (TestParameter)
        tests = {'AcademicTest0','AcademicTest1', 'AcademicTest2', 'AcademicTest3'...
            'AcademicTest4'}
    end

    methods (Test, TestTags = {'Academic'})

        function testFast(testCase, tests)
            s.filename = tests;
            test = AcademicProblem(s);
            test.compute();
            obt = test.result;
            org = load(tests).result;
            err = abs(obt-org);
            tol = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
            close all
        end

    end

end
