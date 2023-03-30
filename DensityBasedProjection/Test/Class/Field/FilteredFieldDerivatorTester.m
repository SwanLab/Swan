classdef FilteredFieldDerivatorTester < handle
    properties (Access = private)
        data
        results
        iterations
        tolerateError
        expectedResult

        filterParameters
    end

    methods (Access = public)
        function obj = FilteredFieldDerivatorTester(iterations)
            obj.iterations = iterations;
            obj.tolerateError = 1e-10;
            obj.loadData()
        end
        function compute(obj)
            % Create the objects

            s.H =obj.filterParameters.H;
            s.Hs =obj.filterParameters.Hs;
            obj.results = FilteredFieldFieldDerivator(s);
            obj.results.compute();
        end
        function loadResults(obj,cParams)
            % In case is testing an external class
            obj.results.derivedFilteredField = cParams.results;
        end 
        function validate(obj)
            % ValidafilterParameterstor
            if abs(obj.expectedResult-obj.results.derivedFilteredField)< obj.tolerateError
                %fprintf('{Stifness matrix}');cprintf('_green', '{true}');disp('|');
                disp('Filtered Field Derivator |OK!|')
            else
                warning('Error in Filtered Field Derivator')
            end
        end
    end
    methods (Access = private)
        function loadData(obj)
            % Load results
            if obj.iterations == 3
                file = fullfile("DensityBasedProjection",'Test','Data','filterParameters');
                s = load(file);
                obj.filterParameters = s.filterParameters;

                file = fullfile("DensityBasedProjection",'Test','Data','derivedFilteredField');
                s = load(file);
                obj.expectedResult = s.derivedFilteredField;

            else
                error('No test Data for the current iterations')
            end
        end
    end
end