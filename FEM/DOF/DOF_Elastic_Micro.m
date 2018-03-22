classdef DOF_Elastic_Micro < DOF_Elastic
    %DOF_Elastic_Micro Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % !! Consider moving periodic to MICRO case !!
        periodic_free % Perioic
        periodic_constrained
    end
    
    methods
        function obj = DOF_Elastic_Micro(problemID,geometry,mesh)
            obj@DOF_Elastic(problemID,geometry,mesh);
            obj.periodic_free = obj.compute_periodic_nodes(obj.master_slave(:,1),obj.nunkn);
            obj.periodic_constrained = obj.compute_periodic_nodes(obj.master_slave(:,2),obj.nunkn);
            obj.constrained{1}=obj.compute_constrained_dof(1);
            obj.free{1}=obj.compute_free_dof(1);
        end
        
        function constrained = compute_constrained_dof(obj,ifield)
            constrained = [obj.periodic_constrained;obj.dirichlet{ifield}];
        end
    end
end

