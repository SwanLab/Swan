classdef AcademicTestsSuite < handle

    methods

        function obj = AcademicTestsSuite()
            path = './Test/AcademicTests/AcademicTests.m';
            suite = matlab.unittest.TestSuite.fromFile(path, 'Tag','Academic');
            results = suite.run;
            table(results)
        end

    end
end