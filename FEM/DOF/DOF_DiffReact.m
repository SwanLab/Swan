classdef DOF_DiffReact < DOF
    %DOF_DiffReact Summary of this class goes here
    %   Detailed explanation goes here
    properties
        master_slave
    end
    
    methods
        function obj = DOF_DiffReact(filename,geometry) % Replace mesh for pdim
            obj.nunkn = 1;
            [dirichlet_data,neumann_data,full_dirichlet_data,master_slave] = Preprocess.getBC_mechanics(filename);
            obj.computeDOF(geometry,dirichlet_data,neumann_data,full_dirichlet_data);
            obj.master_slave = master_slave;
        end
    end
end

