classdef PerformanceTests < handle & matlab.unittest.TestCase

    properties (TestParameter)
        
    end


    methods (Test, TestTags = {'Performance'})

        function test2D(testCase)
            index = 1;
            s.dim    = '2D';
            s.length    = 1;
            s.height = 0.1;
            steps = 0.0001:0.0001:0.1;
            index = 1;
            tic
            % tic should start AFTer assembling the connec
            % maybe store it?
            temps = zeros(size(steps));
            nelems = zeros(size(steps));
%             for step = steps
            step = 0.0001;
                test = PerformanceTest(s);
                sol = test.compute(step);
                nelem = sol.getMesh().nelem;
                disp(nelem)
                temps(index) = toc;
                nelems(index) = nelem;
                index = index+1;
%             end
            disp(temps)
            disp(nelems)
            plot(nelems, temps)
        end

%         function test3D(testCase)
%             index = 1;
%             for i = 0.2:0.1:0.5
%                 gg = @() testCase.example3D(i);
%                 temps(index) = timeit(gg);
%                 connecs(index) = size(testCase.example3D(i),1);
%                 index = index+1;
%             end
%         end

    end

    methods (Access = private)

        function connec = example2D(testCase, step)
            s.dim    = '2D';
            s.length    = 1;
            s.height = 0.1;
            test = PerformanceTest(s);
            sol = test.compute(step);
            connec = test.dofConnec;
        end

        function connec = example3D(testCase, step)
            s.dim    = '3D';
            s.length    = 1;
            s.height = 0.1;
            test = PerformanceTest(s);
            sol = test.compute(step);
            connec = test.dofConnec;
        end

    end

end