classdef CohesiveComputer < handle
    
    properties (Access = public)
        mesh
        boundaryConditions       
        functional
        tolerance
        maxIter
    end
    
    properties (Access = private)
        optimizer
        updater
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = CohesiveComputer(cParams)
            obj.init(cParams)
            obj.setOptimizer(cParams)
        end

        function outputData = compute(obj)
            u = LagrangianFunction.create(obj.mesh,2,'P1');
            cost = 0;

            nSteps = length(obj.boundaryConditions.bcValues);
            for iStep = 1:nSteps
                [u,bc] = obj.preprocess(iStep,nSteps,u);
                [u,F,cost,iterMax] = obj.updater.update(u,bc,cost);
                obj.postprocess(iStep,u,F,cost,iterMax)
            end
            outputData = obj.monitor.data;
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh               = cParams.mesh;
            obj.boundaryConditions = cParams.boundaryConditions;
            obj.functional         = cParams.functional;         
        end
        
        function setOptimizer(obj,cParams)
            s.functional = obj.functional;
            s.tolerance  = cParams.tolerance;
            s.maxIter    = cParams.maxIter;
            s.solverType = cParams.solverType;
            s.monitor    = cParams.monitor;
            obj.updater = DisplacementUpdater(s);
        end
        

        function [u,bc] = preprocess(obj,iStep,nSteps,u)
            % obj.monitor.printStep(iStep,nSteps)
            bc = obj.boundaryConditions.nextStep();
            % u  = obj.computeInitialDisplacement(u,bc);
        end

        


    end
    
end