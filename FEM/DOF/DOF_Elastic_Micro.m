classdef DOF_Elastic_Micro < DOF
    
    
    
    properties
        % !! Consider moving periodic to MICRO case !!
        periodic_free % Perioic
        periodic_constrained
    end
    
    properties (Access = private)
       interpolation
       nFields
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
            obj.nFields = nFields;
            obj.interpolation = interp;            
            [dirichlet_data,neumann_data,full_dirichlet_data,master_slave] = Preprocess.getBC_mechanics(filename);
            if isempty(master_slave)
                mesh.computeMasterSlaveNodes;
                master_slave = mesh.masterSlaveNodes;
            end
            obj.computeData(dirichlet_data,neumann_data,full_dirichlet_data,master_slave,mesh)
        end
        
        function constrained = compute_constrained_dof(obj,ifield)
            constrained = [obj.periodic_constrained;obj.dirichlet{ifield}];
        end
        
        function setMesh(obj,mesh)
            cD = CellNodesDescriptor(mesh.coord);
            corn = cD.cornerNodes;
            
            ncorners = size(corn,1);
            repCorn = repmat(corn,1,2)';
            repDir  = repmat((1:mesh.ndim)',ncorners,1);
            repCorn = repmat(corn,1,2)';
            corners = [repCorn(:), repDir(:), zeros(length(repDir(:)),1)];
            
            msRelator = MasterSlaveRelator(mesh.coord);
            master_slave = msRelator.getRelation();
            dirichlet_data{1} = corners;
            neumann_data = [];
            full_dirichlet_data = [];
            
            obj.computeData(dirichlet_data,neumann_data,full_dirichlet_data,master_slave,mesh)

            
        end
    end
    
    methods (Access = private)
        
        function computeData(obj,dirichlet_data,neumann_data,full_dirichlet_data,master_slave,mesh)
            obj.getDOFconditions(obj.nFields,dirichlet_data,neumann_data,full_dirichlet_data);
            obj.master_slave = master_slave;
            obj.periodic_free = obj.compute_periodic_nodes(obj.master_slave(:,1),obj.nunkn);
            obj.periodic_constrained = obj.compute_periodic_nodes(obj.master_slave(:,2),obj.nunkn);
            obj.computeDOF(mesh,obj.interpolation);
        end
    end
    
    
end

