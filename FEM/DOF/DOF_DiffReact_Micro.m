classdef DOF_DiffReact_Micro < DOF_DiffReact
    %DOF_DiffReact_Micro Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        master_slave
        periodic_free
        periodic_constrained
    end
    
    methods
        function obj = DOF_DiffReact_Micro(problemID,geometry)
            obj@DOF_DiffReact(geometry);
            [~,~,~,obj.master_slave] = Preprocess.getBC_mechanics(problemID);
            obj.periodic_free = obj.compute_periodic_nodes(obj.master_slave(:,1),obj.nunkn);
            obj.periodic_constrained = obj.compute_periodic_nodes(obj.master_slave(:,2),obj.nunkn);
            obj.constrained{1}=obj.compute_constrained_dof(1);
            obj.free{1}=obj.compute_free_dof(1);
        end
        
        function constrained = compute_constrained_dof(obj,ifield)
            % !! Cutrillu !!
            if ~isempty(obj.dirichlet)
                constrained = [obj.periodic_constrained;obj.dirichlet{ifield}];
            else
                constrained = [obj.periodic_constrained;obj.dirichlet];
            end
        end
    end
end

