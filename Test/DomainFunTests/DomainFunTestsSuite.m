classdef DomainFunTestsSuite < handle

    methods

        function obj = DomainFunTestsSuite()
            path = './Test/DomainFunTests/DomainFunTests.m';
            suite = matlab.unittest.TestSuite.fromFile(path, 'Tag','DomainFun');
            results = suite.run;
            table(results)
        end

    end

end