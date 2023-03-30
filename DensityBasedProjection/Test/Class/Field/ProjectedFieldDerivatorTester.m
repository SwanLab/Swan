classdef ProjectedFieldDerivatorTester < handle
    properties (Access = private)
        data
        results
        iterations
        tolerateError

        projectorParameters
        filteredField
        expectedResults
    end

    methods (Access = public)
        function obj = ProjectedFieldDerivatorTester(iterations)
            obj.iterations = iterations;
            obj.tolerateError = 1e-10;
            obj.loadData()
        end
        function compute(obj)
            %% Create the objects
            s.beta =obj.projectorParameters.beta;
            s.filteredField = obj.filteredField;
            s.eta = obj.projectorParameters.eta;

            obj.results = ProjectedFieldFilteredFieldDerivator(s);
            obj.results.compute();

        end
        function loadResults(obj,cParams)
            % In case is testing an external class
            obj.results.derivedProjectedField = cParams.results;
        end 
        function validate(obj)
            %% ValidafilterParameterstor
            if abs(obj.expectedResults-obj.results.derivedProjectedField)< obj.tolerateError
                %fprintf('{Stifness matrix}');cprintf('_green', '{true}');disp('|');
                disp('Projected Field Derivator |OK!|')
            else
                warning('Error in Projected Field Derivator')
            end
        end

    end
    methods (Access = private)
        function loadData(obj)
            %% Load results
            if obj.iterations == 3
                s = fullfile("DensityBasedProjection",'Test','Data','projectorParameters.mat');
                a = load(s);
                obj.projectorParameters = a.projectorParameters;
                
                s = fullfile("DensityBasedProjection",'Test','Data','filteredField.mat');
                a = load(s);
                obj.filteredField = a.filteredField;
                
                s = fullfile("DensityBasedProjection",'Test','Data','derivedProjectedField.mat');
                a = load(s);
                obj.expectedResults = a.derivedProjectedField;
            else
                error('No test Data for the current iterations')
            end
        end

    end
end