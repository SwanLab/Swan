classdef NewTestSuite < handle
    
    properties
        Property1
    end
    
    methods

        function obj = NewTestSuite()
            suite = matlab.unittest.TestSuite.fromFile('NewFemTests.m');
            results = suite.run;
            table(results)
        end

    end
end

