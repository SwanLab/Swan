classdef OptimizerTester < handle
    %UNTITLED40 Summary of this class goes here
    %   Detailed explanation goes here

    properties (Access = private)
        results
        tolerateError
        expectedResult

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
                disp('Step 1: Cost |OK!|')
            else
                warning('Step 1: Error in Cosgt')
            end
            if abs(obj.expectedResult.derivedCost.E-obj.results.E.designCost.derivedCost)< obj.tolerateError
                %fprintf('{Stifness matrix}');cprintf('_green', '{true}');disp('|');
                disp('Step 1: Cost |OK!|')
            else
                warning('Step 1: Error in Cosgt')
            end

            if abs(obj.expectedResult.cost.E-obj.results.E.designCost.cost)< obj.tolerateError
                %fprintf('{Stifness matrix}');cprintf('_green', '{true}');disp('|');
                disp('Step 1: Cost |OK!|')
            else
                warning('Step 1: Error in Cosgt')
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
            
            obj.results = Optimizer(s);
            obj.results.compute

        end
        function loadExpectedResults(obj)
            % Load results
            file = fullfile("DensityBasedProjection",'Test','Data','Optimizer','cost');
            s = load(file);
            obj.expectedResult.cost = s.cost;

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
        end
    end
end