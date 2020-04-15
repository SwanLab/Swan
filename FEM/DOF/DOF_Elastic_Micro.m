classdef DOF_Elastic_Micro < DOF
    
    properties
        % !! Consider moving periodic to MICRO case !!
        periodic_free % Perioic
        periodic_constrained
        master_slave
    end
    
    methods
        function obj = DOF_Elastic_Micro(filename,mesh,pdim,nFields,interp)
            switch pdim
                case '2D'
                    obj.nunkn = 2;
                case '3D'
                    obj.nunkn = 3;
            end
            [dirichlet_data,neumann_data,full_dirichlet_data,master_slave] = Preprocess.getBC_mechanics(filename);
            obj.getDOFconditions(nFields,dirichlet_data,neumann_data,full_dirichlet_data);
            obj.master_slave = master_slave;
            obj.periodic_free = obj.compute_periodic_nodes(obj.master_slave(:,1),obj.nunkn);
            obj.periodic_constrained = obj.compute_periodic_nodes(obj.master_slave(:,2),obj.nunkn);
            obj.computeDOF(mesh,interp);
        end       
        
        function constrained = compute_constrained_dof(obj,ifield)
            constrained = [obj.periodic_constrained;obj.dirichlet{ifield}];
        end
    end
end

