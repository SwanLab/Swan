classdef FieldProjectorTester < handle
    properties (Access = private)
        data
        results
        iterations
        tolerateError
        expectedResults

        projectorParameters
        filteredField

    end

    methods (Access = public)
        function obj = FieldProjectorTester(iterations)
            obj.iterations = iterations;
            obj.tolerateError = 1e-10;
            obj.loadData()
        end
        function compute(obj)
            % Create the objects
            s.beta =obj.projectorParameters.beta;
            s.filteredField = obj.filteredField;
            s.eta = obj.projectorParameters.eta;

            obj.results = FieldProjector(s);
            obj.results.compute();
        end
        function loadResults(obj,cParams)
            % In case is testing an external class
            obj.results.projectedField = cParams.results;
        end 
        function validate(obj)
            % ValidafilterParameterstor
            if abs(obj.expectedResults-obj.results.projectedField)< obj.tolerateError
                %fprintf('{Stifness matrix}');cprintf('_green', '{true}');disp('|');
                disp('Projector |OK!|')
            else
                warning('Error in Projector')
            end
        end
    end
    methods (Access = private)
        function loadData(obj)
            % Load results
            if obj.iterations == 3
                s = fullfile("DensityBasedProjection",'Test','Data','projectorParameters.mat');
                a = load(s);
                obj.projectorParameters = a.projectorParameters;
                
                s = fullfile("DensityBasedProjection",'Test','Data','filteredField.mat');
                a = load(s);
                obj.filteredField = a.filteredField;
                
                s = fullfile("DensityBasedProjection",'Test','Data','projectedField.mat');
                a = load(s);
                obj.expectedResults = a.projectedField;

            else
                error('No test Data for the current iterations')
            end
        end

    end
end