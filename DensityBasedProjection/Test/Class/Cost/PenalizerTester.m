classdef PenalizerTester < handle
    properties (Access = private)
        data
        results
        iterations
        tolerateError
        expectedResult

        projectedField
        structure
    end

    methods (Access = public)
        function obj = PenalizerTester(iterations)
            obj.iterations = iterations;
            obj.tolerateError = 1e-10;
            obj.loadData()
        end
        function compute(obj)
            % Create the objects

            s.elasticModuleMinimun = obj.structure.elasticModuleMinimun;
            s.elasticModuleNeutral = obj.structure.elasticModuleNeutral;
            s.projectedField = obj.projectedField;
            s.penalization = obj.structure.penalization;
            s.nonPenalizedVariable = obj.structure.elementalStiffnessMatrix;
            obj.results  = Penalizer(s);
            obj.results.penalize();
        end
        function loadResults(obj,cParams)
            % In case is testing an external class
            obj.results.filteredField = cParams.results;
        end 
        function validate(obj)
            % ValidafilterParameterstor
            if abs(obj.expectedResult-obj.results.penalizedVariable)< obj.tolerateError
                %fprintf('{Stifness matrix}');cprintf('_green', '{true}');disp('|');
                disp('Penalizer |OK!|')
            else
                warning('Error in Penalizer')
            end
        end
    end
    methods (Access = private)
        function loadData(obj)
            % Load results
            if obj.iterations == 3

                file = fullfile("DensityBasedProjection",'Test','Data','structure');
                s = load(file);
                obj.structure = s.structure;
                file = fullfile("DensityBasedProjection",'Test','Data','projectedField');
                s = load(file);
                obj.projectedField = s.projectedField;
                file = fullfile("DensityBasedProjection",'Test','Data','penalizedVariable');
                s = load(file);
                obj.expectedResult = s.penalizedVariable;

            else
                error('No test Data for the current iterations')
            end
        end
    end
end