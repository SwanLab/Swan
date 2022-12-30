classdef DehomogenizationTestsSuite < handle
    
    methods

        function obj = DehomogenizationTestsSuite()
            warning('off', 'MATLAB:structOnObject')
            results = runtests("DehomogenizingTests","Tag","PlottingTests", 'Debug', true);
            table(results)
            warning('on', 'MATLAB:structOnObject')
        end

    end
end
