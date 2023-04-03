classdef FieldFilterTester < handle
    properties (Access = private)
        data
        results
        iterations
        tolerateError
        expectedResult

        filterParameters
        field
    end

    methods (Access = public)
        function obj = FieldFilterTester(iterations)
            obj.iterations = iterations;
            obj.tolerateError = 1e-10;
            obj.loadData()
        end
        function compute(obj)
            % Create the objects
            s.filterParameters =obj.filterParameters;
            s.field = obj.field;
            obj.results = FieldFilter(s);
            obj.results.compute();
        end
        function loadResults(obj,cParams)
            % In case is testing an external class
            obj.results.filteredField = cParams.results;
        end 
        function validate(obj)
            % ValidafilterParameterstor
            if abs(obj.expectedResult-obj.results.filteredField)< obj.tolerateError
                %fprintf('{Stifness matrix}');cprintf('_green', '{true}');disp('|');
                disp('Filter |OK!|')
            else
                warning('Error in Filter')
            end
        end
    end
    methods (Access = private)
        function loadData(obj)
            % Load results
            if obj.iterations == 3
                file = fullfile("DensityBasedProjection",'Test','Data','filteredField');
                s = load(file);
                obj.expectedResult = s.filteredField;
            else
                error('No test Data for the current iterations')
            end
        end
    end
end