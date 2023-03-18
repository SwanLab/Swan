classdef fieldTester < handle
    properties (Access = private)
        expectedResults
        results
        test
        iterations
        tolerateError
    end

    methods (Access = public)
        function obj = fieldTester(iterations,tolerateError)
            obj.iterations = iterations;
            obj.tolerateError = tolerateError;
            obj.computeTest()
        end
    end
    methods (Access = private)
        function computeTest(obj)
            obj.loadData()
            obj.createComplianceObjects()
            obj.validate()
        end
        function loadData(obj)
            %% Load results
            if obj.iterations == 3
                file = fullfile("DensityBasedProjection",'Test','Data','ResultsData3Iterations.mat');
                obj.expectedResults = load(file);
            elseif obj.iterations == 5
                file = fullfile("DensityBasedProjection",'Test','Data','ResultsData5Iterations.mat');
                obj.expectedResults = load(file);
            else
                error('No test Data for the current iterations')
            end
        end
        function createComplianceObjects(obj)
            %% Create the objects
            s.iterations = obj.iterations;
            obj.results = ComplianceRobustComputer(s);
            obj.results.compute();
        end
        function validate(obj)
            %% Validator
            if abs(obj.expectedResults.results.projectedField.E-obj.results.projectedField.E)< obj.tolerateError
                %fprintf('{Stifness matrix}');cprintf('_green', '{true}');disp('|');
                disp('Projected Field E |OK!|')
            else
                warning('Error in Projected Field E')
            end

            if abs(obj.expectedResults.results.projectedField.I - obj.results.projectedField.I)< obj.tolerateError
                %fprintf('{Stifness matrix}');cprintf('_green', '{true}');disp('|');
                disp('Projected Field I |OK!|')
            else
                warning('Error in Projected Field I')
            end
            if abs(obj.expectedResults.results.projectedField.D - obj.results.projectedField.D)< obj.tolerateError
                %fprintf('{Stifness matrix}');cprintf('_green', '{true}');disp('|');
                disp('Projected Field D |OK!|')
            else
                warning('Error in Projected Field D')
            end
            close all
        end
    end
end