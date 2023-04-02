classdef OptimizerTester < handle
    %UNTITLED40 Summary of this class goes here
    %   Detailed explanation goes here

    properties (Access = private)
        results
        tolerateError
        expectedResult
        data

        E
        I
        D
    end

    methods (Access = public)
        function obj = OptimizerTester()
            obj.tolerateError = 1e-10;
            obj.createObject();
            obj.loadExpectedResults();
        end
        function validate(obj)
            % ValidafilterParameterstor
            if abs(obj.expectedResult.cost.E-obj.results.E.designCost.cost)< obj.tolerateError
                %fprintf('{Stifness matrix}');cprintf('_green', '{true}');disp('|');
                disp('Step 1: cost |OK!|')
            else
                warning('Step 1: Error in Cost')
            end
            if abs(obj.expectedResult.derivedCost.E-obj.results.E.designCost.derivedCost)< obj.tolerateError
                %fprintf('{Stifness matrix}');cprintf('_green', '{true}');disp('|');
                disp('Step 2: derived cost |OK!|')
            else
                warning('Step 2: Error in derived Cost')
            end

            if abs(obj.expectedResult.volumen.volfracD-obj.results.D.designVolumen.volumenFraction)> obj.tolerateError
                %fprintf('{Stifness matrix}');cprintf('_green', '{true}');disp('|');
                warning('Step 3: Error in volumen fraction')
            elseif abs(obj.expectedResult.volumen.value-obj.results.D.designVolumen.volumen)> obj.tolerateError
                warning('Step 3: Error in volumen value')
            elseif abs(obj.expectedResult.volumen.derivated-obj.results.D.designVolumen.derivedVolumen)> obj.tolerateError
                warning('Step 3: Error in derived volumen ')
            else
                disp('Step 3: volumen values |OK!|')                
            end           
            if abs(obj.expectedResult.field-obj.results.E.designField.field)< obj.tolerateError
                %fprintf('{Stifness matrix}');cprintf('_green', '{true}');disp('|');
                disp('Step 4: solver (new field) |OK!|')
            else
                warning('Step 4: Error in solver (new field)')
            end 
            if abs(obj.expectedResult.filteredField-obj.results.E.designField.filteredField)< obj.tolerateError
                %fprintf('{Stifness matrix}');cprintf('_green', '{true}');disp('|');
                disp('Step 5: filtered new field |OK!|')
            else
                warning('Step 5: Error in filtered new field')
            end 
            if abs(obj.expectedResult.projectedField.E-obj.results.E.designField.projectedField)< obj.tolerateError
                %fprintf('{Stifness matrix}');cprintf('_green', '{true}');disp('|');
                disp('Step 6: projected new field |OK!|')
            else
                warning('Step 6: Error in projected new field')
            end
            file = fullfile("DensityBasedProjection",'Test','Data','Optimizer','E');
            s = load(file);
            if abs(s.E.designField.derivedProjectedField-obj.expectedResult.derivedProjectedField.E)< obj.tolerateError
                %fprintf('{Stifness matrix}');cprintf('_green', '{true}');disp('|');
                disp('Step 7: derived projected new field |OK!|')
            else
                warning('Step 7: Error in derived projected new field')
            end    
        end 
    end
    methods (Access = private)
        function createObject(obj)
            file = fullfile("DensityBasedProjection",'Test','Data','Optimizer','E');
            a = load(file);
            s.E = a.E;

            file = fullfile("DensityBasedProjection",'Test','Data','Optimizer','I');
            a = load(file);
            s.I = a.I;

            file = fullfile("DensityBasedProjection",'Test','Data','Optimizer','D');
            a = load(file);
            s.D = a.D;

            file = fullfile("DensityBasedProjection",'Test','Data','mesh');
            a = load(file);
            s.mesh = a.mesh;

            file = fullfile("DensityBasedProjection",'Test','Data','structure');
            a = load(file);
            s.structure = a.structure;
            s.structure.elementType = 'Square';

            file = fullfile("DensityBasedProjection",'Test','Data','projectorParameters');
            a = load(file);
            s.projector = a.projectorParameters;

            file = fullfile("DensityBasedProjection",'Test','Data','filterParameters');
            a = load(file);
            s.filterParameters = a.filterParameters;

            file = fullfile("DensityBasedProjection",'Test','Data','solverParameters');
            a = load(file);
            s.solverParameters = a.solverParameters;

            s.iterations = 1;
            
            obj.data = s;
            obj.results = Optimizer(s);
            obj.results.compute
            
             

        end
        function loadExpectedResults(obj)
            % Load results
            file = fullfile("DensityBasedProjection",'Test','Data','Optimizer','cost');
            s = load(file);
            obj.expectedResult.cost = s.cost;

            file = fullfile("DensityBasedProjection",'Test','Data','Optimizer','field');
            s = load(file);
            obj.expectedResult.field = s.field;

            
            file = fullfile("DensityBasedProjection",'Test','Data','Optimizer','projectedField');
            s = load(file);
            obj.expectedResult.projectedField = s.projectedField;

            file = fullfile("DensityBasedProjection",'Test','Data','Optimizer','filteredField');
            s = load(file);
            obj.expectedResult.filteredField = s.filteredField;

            file = fullfile("DensityBasedProjection",'Test','Data','Optimizer','costChange');
            s = load(file);
            obj.expectedResult.costChange = s.costChange;

            file = fullfile("DensityBasedProjection",'Test','Data','Optimizer','derivedCost');
            s = load(file);
            obj.expectedResult.derivedCost = s.derivedCost;

            file = fullfile("DensityBasedProjection",'Test','Data','Optimizer','derivedProjectedField');
            s = load(file);
            obj.expectedResult.derivedProjectedField = s.derivedProjectedField;

            file = fullfile("DensityBasedProjection",'Test','Data','Optimizer','displacement');
            s = load(file);
            obj.expectedResult.displacement = s.displacement;

            file = fullfile("DensityBasedProjection",'Test','Data','Optimizer','volumen');
            s = load(file);
            obj.expectedResult.volumen = s.volumen;


        end
    end
end