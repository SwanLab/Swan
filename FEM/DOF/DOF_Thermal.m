classdef DOF_Thermal < DOF
    %DOF_Thermal Summary of this class goes here
    %   Detailed explanation goes here
    
    methods
        function obj = DOF_Thermal(filename,geometry) % Replace mesh for pdim
            obj.nunkn = 1;
            [dirichlet_data,neumann_data,full_dirichlet_data,master_slave] = Preprocess.getBC_mechanics(filename);
            obj.computeDOF(geometry,dirichlet_data,neumann_data,full_dirichlet_data,master_slave);
        end
    end
end

