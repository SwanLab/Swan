classdef MicroShapeTestsSuite < handle

    methods

        function obj = MicroShapeTestsSuite()
            path = './Test/MicroShapeTests/MicroShapeTests.m';
            suite = matlab.unittest.TestSuite.fromFile(path, 'Tag','MicroShape');
            results = suite.run;
            table(results)
        end

    end
end