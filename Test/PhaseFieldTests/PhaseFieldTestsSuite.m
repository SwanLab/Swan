classdef PhaseFieldTestsSuite < handle

    methods

        function obj = PhaseFieldTestsSuite()
            path = './Test/PhaseFieldTests/PhaseFieldTests.m';
            suite = matlab.unittest.TestSuite.fromFile(path, 'Tag','PF');
            results = suite.run;
            table(results)
        end

    end

end