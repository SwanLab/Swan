classdef CostComputerTester < handle
    properties (Access = private)
        iterations
        tolerateError
        cost

        designCost        
    end

    methods (Access = public)
        function obj = CostComputerTester()
            obj.iterations = 3;
            obj.tolerateError = 1e-10;
            obj.createDesignCost();
        end
        function testFEM(obj)
            obj.designCost.computeCost();
            results = obj.designCost.cost;
            %Validate Results
            s.results = results;
            B = FEMcomputerTester(3);
            B.loadResults(s);
            B.validate;
        end
        function testCostDerivator(obj)
            obj.testFEM()
            obj.designCost.deriveCost();
            results = obj.designCost.derivedCost;
            %Validate Results
            s.results = results;
            B = CostFieldDerivatorTester(3);
            B.loadResults(s);
            B.validate;
        end

    end
    methods (Access = private)
        function createDesignCost(obj)
            %Load Initial Paramaters

            file = fullfile("DensityBasedProjection",'Test','Data','structure.mat');            
            a = load(file);
            s.structure = a.structure;

            file = fullfile("DensityBasedProjection",'Test','Data','mesh.mat');            
            a = load(file);
            s.mesh = a.mesh;

            file = fullfile("DensityBasedProjection",'Test','Data','designFields.mat');            
            a = load(file);
            s.designFields = a.designFields;

            file = fullfile("DensityBasedProjection",'Test','Data','filterParameters.mat');            
            a = load(file);
            s.filterParameters = a.filterParameters;

            %Create the DesignFieldObject
            obj.designCost = CostComputer(s);
        end 

    end
end
