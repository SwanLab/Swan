classdef OptimizerTester < handle
    properties (Access = private)
        data
        results
        iterations
        tolerateError
        expectedResult

        filterParameters
        projectParameters
        
        field
        mesh
        volumen 
        solverParameters
        cost
        structure
        filteredField
        projectedField
    end

    methods (Access = public)
        function obj = OptimizerTester(iterations)
            obj.iterations = iterations;
            obj.tolerateError = 1e-10;
            obj.loadData()
        end
        function compute(obj)
            % Create the objects

            s.mesh = obj.mesh;
            s.structure=obj.structure;
            s.structure.elementType = 'Square';
            s.projector = obj.projectParameters;
            s.filterParameters = obj.filterParameters;
            s.cost = obj.cost;
            s.cost.initial = obj.cost.E;
            s.solverParameters =obj.solverParameters;
            s.iterations = obj.iterations;
            s.field = obj.field;
            s.filteredField = obj.filteredField;
            s.projectedField.E = obj.projectedField;
            s.volumen = obj.volumen;
            obj.results = Optimizer(s);
            obj.results.compute();
        end
        function loadResults(obj,cParams)
            % In case is testing an external class
            obj.results.filteredField = cParams.results;
        end 
        function validate(obj)
            % ValidafilterParameterstor
            if abs(obj.expectedResult-obj.results.filteredField)< obj.tolerateError
                %fprintf('{Stifness matrix}');cprintf('_green', '{true}');disp('|');
                disp('Filter |OK!|')
            else
                warning('Error in Filter')
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
                file = fullfile("DensityBasedProjection",'Test','Data','structure');
                s = load(file);
                obj.structure = s.structure;                      
                file = fullfile("DensityBasedProjection",'Test','Data','filterParameters');
                s = load(file);
                obj.filterParameters = s.filterParameters;
                file = fullfile("DensityBasedProjection",'Test','Data','projectorParameters');
                s = load(file);
                obj.projectParameters = s.projectorParameters;  
                file = fullfile("DensityBasedProjection",'Test','Data','cost');
                s = load(file);
                obj.cost = s.cost;
                file = fullfile("DensityBasedProjection",'Test','Data','field');
                s = load(file);
                obj.field = s.field;
                file = fullfile("DensityBasedProjection",'Test','Data','filteredField');
                s = load(file);
                obj.filteredField = s.filteredField;
                file = fullfile("DensityBasedProjection",'Test','Data','solverParameters');
                s = load(file);                
                obj.solverParameters = s.solverParameters;
                file = fullfile("DensityBasedProjection",'Test','Data','projectedField');
                s = load(file);                
                obj.projectedField = s.projectedField;                
                file = fullfile("DensityBasedProjection",'Test','Data','volumen');
                s = load(file);
                obj.volumen = s.volumen;                
            else
                error('No test Data for the current iterations')
            end
        end
    end
end

