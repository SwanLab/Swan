classdef PhaseFieldTestsSuite < handle

    methods

        function obj = PhaseFieldTestsSuite()
            path = './tests/Source/PhaseFieldTests/PhaseFieldTests.m';
            suite = matlab.unittest.TestSuite.fromFile(path,'Tag','PF');
            results = suite.run;
            table(results)
        end

    end

end