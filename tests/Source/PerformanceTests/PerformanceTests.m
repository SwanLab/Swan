classdef PerformanceTests < handle & matlab.unittest.TestCase

    properties (TestParameter)
        
    end


    methods (Test, TestTags = {'Performance'})

        function test2D(testCase)
            index = 1;
            for i = 0.01:0.005:0.5
                gg = @() testCase.example2D(i);
                temps(index) = timeit(gg);
                connecs(index) = size(testCase.example2D(i),1);
                index = index+1;
            end
        end

        function test3D(testCase)
            index = 1;
            for i = 0.2:0.1:0.5
                gg = @() testCase.example3D(i);
                temps(index) = timeit(gg);
                connecs(index) = size(testCase.example3D(i),1);
                index = index+1;
            end
        end

    end

    methods (Access = private)

        function connec = example2D(testCase, step)
            s.dim    = '2D';
            s.length    = 1;
            s.height = 0.1;
            test = PerformanceTest(s);
            sol = test.compute(step);
            connec = test.getConnec();
        end

        function connec = example3D(testCase, step)
            s.dim    = '3D';
            s.length    = 1;
            s.height = 0.1;
            test = PerformanceTest(s);
            sol = test.compute(step);
            connec = test.getConnec();
        end

    end

end