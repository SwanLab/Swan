classdef testFemComputation < handle
    
    properties (Abstract, Access = protected)
        testName
        computedVar
    end
    
    properties (Access = protected)
       fem
    end
    
    methods (Access = protected)
        
        function obj = testFemComputation()
           obj.computeVariableThroughFemSolver()
           obj.selectComputedVar();
        end
    end
    
    methods (Access = protected)
        
        function computeVariableThroughFemSolver(obj)
            femSolver = FEM.create(obj.testName);
            femSolver.preProcess;
            femSolver.computeVariables;
            obj.fem = femSolver;
        end
        
    end
    
    methods (Abstract, Access = protected)
        selectComputedVar(obj)
    end
    
end

