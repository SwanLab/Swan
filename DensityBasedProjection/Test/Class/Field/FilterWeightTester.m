classdef FilterWeightTester < handle
    properties (Access = private)
        data
        results
        iterations
        tolerateError
    end

    methods (Access = public)
        function obj = FilterWeightTester(iterations)
            obj.iterations = iterations;
            obj.tolerateError = 1e-10;
            obj.loadData()
        end
        function compute(obj)
            %Create the objects
            s.elementNumberX =obj.data.WeightRatioData.elementNumberX;
            s.elementNumberY = obj.data.WeightRatioData.elementNumberY;
            s.minimunInfluenceRadios = obj.data.WeightRatioData.minimunInfluenceRadios;
            
            obj.results = weightFilterComputer(s);
            obj.results.compute();

        end
        function loadResults(obj,cParams)
            % In case is testing an external class
            obj.results.filteredField = cParams.results;
        end 
        function validate(obj)
            % ValidafilterParameterstor
            if abs(obj.data.WeightRatioData.H-obj.results.H)< obj.tolerateError & abs(obj.data.WeightRatioData.Hs-obj.results.Hs)< obj.tolerateError 
                %fprintf('{Stifness matrix}');cprintf('_green', '{true}');disp('|');
                disp('Filter Weight Tester |OK!|')
            else
                warning('Error in Filter Weight')
            end
            close all
        end

    end
    methods (Access = private)
        function loadData(obj)
            % Load results
            if obj.iterations == 3
                file = fullfile("DensityBasedProjection",'Test','Data','data&resultsWeightRatio.mat');
                obj.data = load(file);
            else
                error('No test Data for the current iterations')
            end
        end
    end
end