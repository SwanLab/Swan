classdef ContinuumDamageTestsSuite < handle

    methods

        function obj = ContinuumDamageTestsSuite()
            path = './Test/ContinuumDamageTests/ContinuumDamageTests.m';
            suite = matlab.unittest.TestSuite.fromFile(path, 'Tag','CD');
            results = suite.run;
            table(results)
        end

    end

end