classdef ReadingFilesTestsSuite < handle
    
    methods

        function obj = ReadingFilesTestsSuite()
            path = './tests/Source/ReadingFilesTests/ReadingFilesTests.m';
            suite = matlab.unittest.TestSuite.fromFile(path, 'Tag','ReadingFiles');
            results = suite.run;
            table(results)
        end

    end
end