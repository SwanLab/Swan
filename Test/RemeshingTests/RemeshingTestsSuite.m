classdef RemeshingTestsSuite < handle
    
    methods

        function obj = RemeshingTestsSuite()
            warning('off', 'MATLAB:structOnObject')
            results = runtests("RemeshingTests","Tag","Remesh", 'Debug', false);
            table(results)
            warning('on', 'MATLAB:structOnObject')
        end

    end
end
