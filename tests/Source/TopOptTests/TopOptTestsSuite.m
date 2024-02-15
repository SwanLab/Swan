classdef TopOptTestsSuite < handle
    
    methods

        function obj = TopOptTestsSuite()
            warning('off', 'MATLAB:structOnObject')
            % testFastDisplacement, testMacro, testMicro
            results = runtests("TopOptTests","ProcedureName","testFastDisplacement", 'Debug', true);
%             results = runtests("TopOptTests","ProcedureName","testAnalyticVsRegularizedPerimeter", 'Debug', true);
%             results = runtests("TopOptTests", 'Debug', true);
            table(results)
            warning('on', 'MATLAB:structOnObject')
        end

    end
end
