classdef ProjectorsTestsSuite < handle
    
    methods

        function obj = ProjectorsTestsSuite()
            warning('off', 'MATLAB:structOnObject')
            results = runtests("ProjectorsTests","ProcedureName","testFastDisplacement", 'Debug', true);
            table(results)
            warning('on', 'MATLAB:structOnObject')
        end

    end
end
