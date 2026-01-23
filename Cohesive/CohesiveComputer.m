classdef CohesiveComputer < handle
    
    properties (Access = public)
        mesh
        boundaryConditions       
        functional
        tolerance
        maxIter
    end
    
    properties (Access = private)
        monitor
        optimizer
        updater

    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = CohesiveComputer(cParams)
            obj.init(cParams)
            obj.setMonitoring(cParams)
            obj.setOptimizer
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh               = cParams.mesh;
            obj.boundaryConditions = cParams.boundaryConditions;
            obj.internalDamageVariable = cParams.internalDamageVariable;
            obj.functional         = cParams.functional;
            obj.tolerance          = cParams.tolerance;
            obj.maxIter            = cParams.maxIter;            
        end
        
        function setOptimizer(obj,cParams)
            s.functional = obj.functional;
            s.monitor    = obj.monitor;
            s.tolerance  = cParams.tolerance;
            s.maxIter    = cParams.maxIter;
            s.solverType = cParams.solverType;
            obj.updater = CohesiveDisplacementUpdater(s);
        end
        
        function setMonitoring(obj,cParams)
            
        end

        
        
        
    end
    
end