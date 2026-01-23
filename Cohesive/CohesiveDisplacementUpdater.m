classdef CohesiveDisplacementUpdater < handle

    
    properties (Access = private)
        functional
            monitor
            mesh
    
            solverType
            solver
    
            tol
            maxIter   
    end

    methods (Access = public)
        
        function obj = CohesiveDisplacementUpdater(cParams)
            obj.init(cParams)
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
                obj.functional = cParams.functional;
                obj.monitor    = cParams.monitor;
                obj.tol        = cParams.tolerance;
                obj.maxIter    = cParams.maxIter;
                obj.solverType = cParams.solverType;
                obj.setSolver();
        end


        function setSolver(obj)
            switch obj.solverType
                case 'Newton'
                    obj.solver = Newton();
                case 'AdaptiveGradient'
                    s.tauMax        = 1e10;
                    s.tau           = 100;
                    s.functional    = obj.functional;
                    obj.solver      = AdaptiveGradient(s);
            end
        end


        function [e,cost] = computeCost(obj,u,bc,costOld)
            cost = obj.functional.computeCost(u,bc);
            e = cost - costOld;
        end










        
    end
    
end