classdef ForceComputerTester < handle
    properties (Access = private)
        data
        results
        iterations
        tolerateError
        expectedResult
        
        mesh

    end

    methods (Access = public)
        function obj = ForceComputerTester(iterations)
            obj.iterations = iterations;
            obj.tolerateError = 1e-10;
            obj.loadData()
        end
        function compute(obj)
            % Create the objects

            s.neumanCondition = obj.mesh.neumanCondition;
            s.output = obj.mesh.output;
            s.elementNumberX = obj.mesh.elementNumberX;
            s.elementNumberY = obj.mesh.elementNumberY;
            obj.results = ForceComputer(s);
            obj.results.compute();
        end
        function loadResults(obj,cParams)
            % In case is testing an external class
            obj.results.filteredField = cParams.results;
        end 
        function validate(obj)
            % ValidafilterParameterstor
            if abs(obj.expectedResult-obj.results.force)< obj.tolerateError
                %fprintf('{Stifness matrix}');cprintf('_green', '{true}');disp('|');
                disp('Force |OK!|')
            else
                warning('Error in Force')
            end
        end
    end
    methods (Access = private)
        function loadData(obj)
            % Load results
            if obj.iterations == 3
                file = fullfile("DensityBasedProjection",'Test','Data','mesh');
                s = load(file);
                obj.mesh = s.mesh;

                file = fullfile("DensityBasedProjection",'Test','Data','force');
                s = load(file);
                obj.expectedResult = s.force;

            else
                error('No test Data for the current iterations')
            end
        end
    end
end