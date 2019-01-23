classdef DOF_Thermal < DOF
    %DOF_Thermal Summary of this class goes here
    %   Detailed explanation goes here
    
    methods
        function obj = DOF_Thermal(filename,geometry) % Replace mesh for pdim
            obj.nunkn = 1;
            [dirichlet_data,neumann_data,full_dirichlet_data] = Preprocess.getBC_mechanics(filename);
            obj.getDOFcondition(geometry,dirichlet_data,neumann_data,full_dirichlet_data);
            obj.computeDOF(geometry);
        end
    end
end

