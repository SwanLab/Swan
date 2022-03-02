classdef TopOptTestsSuite < handle
    
    methods

        function obj = TopOptTestsSuite()
            warning('off', 'MATLAB:structOnObject')
            path = './tests/Source/TopOptTests/TopOptTests.m';
%             suite = matlab.unittest.TestSuite.fromFile(path, 'Tag','FastDisp');%ToPass
%             results = suite.run;
            results = runtests("TopOptTests","ProcedureName","testFastDisplacement", 'Debug', true);
            table(results)
            warning('on', 'MATLAB:structOnObject')
        end

    end
end