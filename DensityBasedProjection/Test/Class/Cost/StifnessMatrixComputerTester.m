classdef StifnessMatrixComputerTester < handle
    properties (Access = private)
        data
        results
        iterations
        tolerateError
        expectedResult
        
        structure
        mesh
        projectedField

    end

    methods (Access = public)
        function obj = StifnessMatrixComputerTester(iterations)
            obj.iterations = iterations;
            obj.tolerateError = 1e-10;
            obj.loadData()
        end
        function compute(obj)
            % Create the objects

            s.structure = obj.structure;
            s.mesh = obj.mesh;
            s.projectedField            =  obj.projectedField;
            obj.results = StifnessMatrixComputer(s);
            obj.results.compute();
        end
        function loadResults(obj,cParams)
            % In case is testing an external class
            obj.results.filteredField = cParams.results;
        end 
        function validate(obj)
            % ValidafilterParameterstor
            if abs(obj.expectedResult-obj.results.globalStifnessMatrix)< obj.tolerateError
                %fprintf('{Stifness matrix}');cprintf('_green', '{true}');disp('|');
                disp('Stifness Matrix Computer |OK!|')
            else
                warning('Error in Stifness Matrix Computer')
            end
        end
    end
    methods (Access = private)
        function loadData(obj)
            % Load results
            if obj.iterations == 3
                file = fullfile("DensityBasedProjection",'Test','Data','structurePenalized');
                s = load(file);
                obj.structure = s.structurePenalized;
                file = fullfile("DensityBasedProjection",'Test','Data','mesh');
                s = load(file);
                obj.mesh = s.mesh;
                file = fullfile("DensityBasedProjection",'Test','Data','projectedField');
                s = load(file);
                obj.projectedField = s.projectedField;

                file = fullfile("DensityBasedProjection",'Test','Data','globalStifnessMatrix');
                s = load(file);
                obj.expectedResult = s.globalStifnessMatrix;


            else
                error('No test Data for the current iterations')
            end
        end
    end
end