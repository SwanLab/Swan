classdef DehomogenizationTestsSuite < handle
    
    methods

        function obj = DehomogenizationTestsSuite()
            warning('off', 'MATLAB:structOnObject')
            results = runtests("DehomogenizingTests","Tag","Singularities", 'Debug', true);
            table(results)
            warning('on', 'MATLAB:structOnObject')
        end

    end
end
