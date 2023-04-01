classdef VolumenComputerTester < handle
    properties (Access = private)
        data
        results
        iterations
        tolerateError
        expectedResult

    end

    methods (Access = public)
        function obj = VolumenComputerTester(iterations)
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
            if abs(obj.expectedResult.volumen-obj.results)< obj.tolerateError
                %fprintf('{Stifness matrix}');cprintf('_green', '{true}');disp('|');
                disp('Volumen Computer |OK!|')
            else
                warning('Error in Volumen Computer')
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