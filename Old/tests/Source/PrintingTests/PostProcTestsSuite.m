classdef PostProcTestsSuite < handle
    
    methods

        function obj = PostProcTestsSuite()
            warning('off', 'MATLAB:structOnObject')
            results = runtests("PostProcTests","ProcedureName","testFastDisplacement", 'Debug', true);
            table(results)
            warning('on', 'MATLAB:structOnObject')
        end

    end
end