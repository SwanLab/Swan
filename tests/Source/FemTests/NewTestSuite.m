classdef NewTestSuite < handle
    
    methods

        function obj = NewTestSuite()
            suite = matlab.unittest.TestSuite.fromFile('NewFemTests.m');
            results = suite.run;
            table(results)
        end

    end
end

