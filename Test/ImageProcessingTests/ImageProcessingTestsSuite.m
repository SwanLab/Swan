classdef ImageProcessingTestsSuite < handle
    
    methods

        function obj = ImageProcessingTestsSuite()
            path = './Test/ImageProcessingTests/ImageProcessingTests.m';
            suite = matlab.unittest.TestSuite.fromFile(path, 'Tag','ImageProcessing');
            results = suite.run;
            table(results)
        end

    end
end