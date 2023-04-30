classdef GlobalTest < handle
    properties (Access = private)
        data
        results
        iterations
        tolerateError
        expectedResult
        Test1

    end

    methods (Access = public)
        function obj = GlobalTest(iterations)
            obj.iterations = iterations;
            obj.tolerateError = 1e-10;
            obj.loadData()
            obj.compute()
            obj.validate()
        end
        
    end
    methods (Access = private)
        function loadData(obj)
            if obj.iterations == 3
                file = fullfile("DensityBasedProjection",'Test','Data','ResultsData3Iterations.mat');
                s = load(file);
                obj.expectedResult = s.results;
            elseif obj.iterations == 5
                file = fullfile("DensityBasedProjection",'Test','Data','ResultsData5Iterations.mat');
                s = load(file);
                obj.expectedResult = s.results;
            else
                error('No test Data for the current iterations')
            end
        end
        function compute(obj)
            % Create the objects
            s.iterations = obj.iterations;
            s.settingsName = 'test_arturo';
            obj.Test1 = ComplianceRobustComputer(s);
            obj.Test1.compute();
        end
        function validate(obj)
            % ValidafilterParameterstor
            if abs(obj.expectedResult.projectedField.E-obj.Test1.E.designField.projectedField)< obj.tolerateError
                %fprintf('{Stifness matrix}');cprintf('_green', '{true}');disp('|');
                disp('Projected Field E |OK!|')
            else
                warning('Error in Projected Field E')
            end

            if abs(obj.expectedResult.projectedField.I - obj.Test1.I.designField.projectedField)< obj.tolerateError
                %fprintf('{Stifness matrix}');cprintf('_green', '{true}');disp('|');
                disp('Projected Field I |OK!|')
            else
                warning('Error in Projected Field I')
            end
            if abs(obj.expectedResult.projectedField.D - obj.Test1.D.designField.projectedField)< obj.tolerateError
                %fprintf('{Stifness matrix}');cprintf('_green', '{true}');disp('|');
                disp('Projected Field D |OK!|')
            else
                warning('Error in Projected Field D')
            end
        end
    end
end

