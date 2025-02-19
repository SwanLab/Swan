classdef TopOptTestsSuite < handle
    
    methods

        function obj = TopOptTestsSuite()
            warning('off', 'MATLAB:structOnObject')
            results = runtests("TopOptTests","ProcedureName","testFastDisplacement", 'Debug', false);
            table(results)
            warning('on', 'MATLAB:structOnObject')
        end

    end
end
