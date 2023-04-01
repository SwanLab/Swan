classdef FEMcomputerTester < handle
    properties (Access = private)
        data
        results
        iterations
        tolerateError
        expectedResult

        projectedField
        structure
        mesh
    end

    methods (Access = public)
        function obj = FEMcomputerTester(iterations)
            obj.iterations = iterations;
            obj.tolerateError = 1e-10;
            obj.loadData()
        end
        function loadResults(obj,cParams)
            % In case is testing an external class
            obj.results.cost = cParams.results;
        end 
        function validate(obj)
            % ValidafilterParameterstor
            if abs(obj.expectedResult.cost-obj.results.cost)< obj.tolerateError
                %fprintf('{Stifness matrix}');cprintf('_green', '{true}');disp('|');
                disp('FEM computer |OK!|')
            else
                warning('Error in FEM computer')
            end
        end
    end
    methods (Access = private)
        function loadData(obj)
            % Load results
            if obj.iterations == 3
                file = fullfile("DensityBasedProjection",'Test','Data','cost');
                s = load(file);
                obj.expectedResult.cost = s.cost.E;               
            else
                error('No test Data for the current iterations')
            end
        end
    end
end