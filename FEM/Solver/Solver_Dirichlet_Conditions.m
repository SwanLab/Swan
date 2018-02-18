classdef Solver_Dirichlet_Conditions < Solver
    %AnalyticalSolver Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        fixnodes
    end
    
    
    methods (Access = public)
        
        function obj = Solver_Dirichlet_Conditions()
        end
        
        
        
        
        function setSolverVariables(obj,data)
            obj.fixnodes = data.fixnodes;
        end
        
    end
    
    
end

