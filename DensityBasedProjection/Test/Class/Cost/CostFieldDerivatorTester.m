classdef CostFieldDerivatorTester < handle
    properties (Access = private)
        data
        results
        iterations
        tolerateError
        expectedResult

        mesh
        structure
        filterParameters
        
        projectorParameters
        cost
        displacement
        
        derivedProjectedField
        designFields

        
    end

    methods (Access = public)
        function obj = CostFieldDerivatorTester(iterations)
            obj.iterations = iterations;
            obj.tolerateError = 1e-10;
            obj.loadData()
        end
        function compute(obj)
            % Create the objects

            s.mesh = obj.mesh;
            s.structure =obj.structure;
            s.filterParameters =obj.filterParameters ;

            s.designFields =obj.designFields ;
            s.projectorParameters=obj.projectorParameters;
            s.cost =obj.cost ;            
            s.displacement =obj.displacement;


            obj.results = CostFieldDerivator(s);
            obj.results.compute();
        end
        function loadResults(obj,cParams)
            % In case is testing an external class
            obj.results.derivedCost = cParams.results;
        end 
        function validate(obj)
            % ValidafilterParameterstor
            if abs(obj.expectedResult-obj.results.derivedCost)< obj.tolerateError
                %fprintf('{Stifness matrix}');cprintf('_green', '{true}');disp('|');
                disp('Cost respected Field Derivator |OK!|')
            else
                warning('Error in Cost respected Field Derivator ')
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
                file = fullfile("DensityBasedProjection",'Test','Data','designFields');
                s = load(file);
                obj.designFields = s.designFields;
                file = fullfile("DensityBasedProjection",'Test','Data','cost');
                s = load(file);
                obj.cost = s.cost;
                           
                file = fullfile("DensityBasedProjection",'Test','Data','cost');
                s = load(file);
                obj.cost = s.cost;
                file = fullfile("DensityBasedProjection",'Test','Data','displacement');
                s = load(file);
                obj.displacement = s.displacement;
                file = fullfile("DensityBasedProjection",'Test','Data','derivedCost');
                s = load(file);
                obj.expectedResult = s.derivedCost;

            else
                error('No test Data for the current iterations')
            end
        end
    end
end