classdef DOF_DiffReact_Micro < DOF_DiffReact
    %DOF_DiffReact_Micro Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % !! Consider moving periodic to MICRO case !!
        periodic_free % Perioic
        periodic_constrained
    end
    
    methods
        function obj = DOF_DiffReact_Micro(problemID,geometry)
            obj@DOF_DiffReact(problemID,geometry);
            obj.periodic_free = obj.compute_periodic_nodes(obj.master_slave(:,1),obj.nunkn);
            obj.periodic_constrained = obj.compute_periodic_nodes(obj.master_slave(:,2),obj.nunkn);
        end
        
        function constrained = compute_constrained_dof(obj,ifield)
            if ~isempty(obj.dirichlet)
                constrained = [obj.periodic_constrained;obj.dirichlet{ifield}];
            else
                constrained = [obj.periodic_constrained;obj.dirichlet];
            end
        end
    end
end

