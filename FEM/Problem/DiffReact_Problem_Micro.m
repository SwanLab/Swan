classdef DiffReact_Problem_Micro < DiffReact_Problem
    %DiffReact_Problem_Micro Summary of this class goes here
    % Detailed explanation goes here
    
    %% Public GetAccess properties definition =============================
    properties (GetAccess = public, SetAccess = public)
    end
    
    %% Restricted properties definition ===================================
    properties %(GetAccess = {?Postprocess,?Physical_Problem_Micro}, SetAccess = protected)
    end
    
    %% Public methods definition ==========================================
    methods (Access = public)
        function obj = DiffReact_Problem_Micro(problemID)
            obj@DiffReact_Problem(problemID);
            obj.dof = DOF_DiffReact_Micro(problemID,obj.geometry);
        end
        
        function preProcess(obj)
            obj.element = Element_DiffReact_Micro(obj.mesh,obj.geometry,obj.material,obj.dof);
            obj.solver = Solver.create;
        end
    end
end

