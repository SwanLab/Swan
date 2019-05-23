classdef testFemComputation < handle
    
    properties (Abstract, Access = protected)
        testName
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
            props.kappa = .9107;
            props.mu    = .3446;        
            femSolver.setMatProps(props);            
            femSolver.computeVariables;
            obj.fem = femSolver;
        end
        
    end
    
    methods (Abstract, Access = protected)
        selectComputedVar(obj)
    end
    
end

