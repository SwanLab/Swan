classdef VolumenFractionTester < handle
    properties (Access = private)
        data
        results
        iterations
        tolerateError
        expectedResult

    end

    methods (Access = public)
        function obj = VolumenFractionTester(iterations)
            obj.iterations = iterations;
            obj.tolerateError = 1e-10;
            obj.loadData()
        end
        function loadResults(obj,cParams)
            % In case is testing an external class
            obj.results = cParams.results;
        end 
        function validate(obj)
            % ValidafilterParameterstor
            if abs(obj.expectedResult.volumenFraction-obj.results)< obj.tolerateError
                %fprintf('{Stifness matrix}');cprintf('_green', '{true}');disp('|');
                disp('Volumen Fraction Computer |OK!|')
            else
                warning('Error in Volumen Fraction Computer')
            end
        end
    end
    methods (Access = private)
        function loadData(obj)
            % Load results
            if obj.iterations == 3
                file = fullfile("DensityBasedProjection",'Test','Data','designVolumen');
                s = load(file);
                obj.expectedResult = s.designVolumen;
            else
                error('No test Data for the current iterations')
            end
        end
    end
end