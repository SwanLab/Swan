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
    methods (Access = protected)
        function setElement(obj)
            obj.element = Element_DiffReact_Micro(obj.mesh,obj.geometry,obj.material,obj.dof);
        end
        
        function setDOFs(obj)
            % !! Coupled to GiD File, pending to set free !!
            obj.dof = DOF_DiffReact_Micro(obj.problemID,obj.geometry);
        end
    end
end
