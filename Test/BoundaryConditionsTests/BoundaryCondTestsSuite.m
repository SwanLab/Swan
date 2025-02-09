classdef BoundaryCondTestsSuite < handle & matlab.unittest.TestCase
    
    methods (Access = public)
        function obj = BoundaryCondTestsSuite()
            path = './Test/BoundaryConditionsTests/BoundaryCondTests.m';
            suite = matlab.unittest.TestSuite.fromFile(path, 'Tag', 'Monolithic' );
            results = suite.run;
            table(results)
        end
    end
end