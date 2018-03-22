classdef DOF_Elastic < DOF
    %DOF_Elastic Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        master_slave
    end
    
    methods
        function obj = DOF_Elastic(filename,geometry,mesh) % Replace mesh for pdim
            switch mesh.pdim
                case '2D'
                    obj.nunkn = 2;
                case '3D'
                    obj.nunkn = 3;
            end
            [dirichlet_data,neumann_data,full_dirichlet_data,master_slave] = Preprocess.getBC_mechanics(filename);
            obj.computeDOF(geometry,dirichlet_data,neumann_data,full_dirichlet_data);
            obj.master_slave = master_slave;
        end
    end
end

