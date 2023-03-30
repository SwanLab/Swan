classdef DisplacementComputerTester < handle
    properties (Access = private)
        data
        results
        iterations
        tolerateError
        expectedResult
        
        force
        mesh
        globalStifnessMatrix

    end

    methods (Access = public)
        function obj = DisplacementComputerTester(iterations)
            obj.iterations = iterations;
            obj.tolerateError = 1e-10;
            obj.loadData()
        end
        function compute(obj)
            % Create the objects

            s.force = obj.force;
            s.globalStifnessMatrix =  obj.globalStifnessMatrix;
            s.freeDegress = obj.mesh.degress.free;
            s.elementNumberX = obj.mesh.elementNumberX;
            s.elementNumberY = obj.mesh.elementNumberY;

            obj.results = DisplacementComputer(s);
            obj.results.compute();
        end
        function loadResults(obj,cParams)
            % In case is testing an external class
            obj.results.filteredField = cParams.results;
        end 
        function validate(obj)
            % ValidafilterParameterstor
            if abs(obj.expectedResult-obj.results.displacement)< obj.tolerateError
                %fprintf('{Stifness matrix}');cprintf('_green', '{true}');disp('|');
                disp('Displacement Computer |OK!|')
            else
                warning('Error in Displacement Computer ')
            end
        end
    end
    methods (Access = private)
        function loadData(obj)
            % Load results
            if obj.iterations == 3
                file = fullfile("DensityBasedProjection",'Test','Data','force');
                s = load(file);
                obj.force = s.force;
                file = fullfile("DensityBasedProjection",'Test','Data','mesh');
                s = load(file);
                obj.mesh = s.mesh;
                file = fullfile("DensityBasedProjection",'Test','Data','globalStifnessMatrix');
                s = load(file);
                obj.globalStifnessMatrix = s.globalStifnessMatrix;

                file = fullfile("DensityBasedProjection",'Test','Data','displacement');
                s = load(file);
                obj.expectedResult = s.displacement;





            else
                error('No test Data for the current iterations')
            end
        end
    end
end