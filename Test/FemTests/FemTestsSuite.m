classdef FemTestsSuite < handle

    methods

        function obj = FemTestsSuite()
            path = './Test/FemTests/FemTests.m';
            suite = matlab.unittest.TestSuite.fromFile(path, 'Tag', 'FEM');
            
            % Create a runner with text output
            runner = matlab.unittest.TestRunner.withTextOutput;
            
            % Add StopOnFailuresPlugin
            import matlab.unittest.plugins.StopOnFailuresPlugin
            runner.addPlugin(StopOnFailuresPlugin);
            
            % Run the suite
            results = runner.run(suite);
            table(results)
        end

    end
end