classdef DOF_Stokes < DOF
    %DOF_Stokes Summary of this class goes here
    %   Detailed explanation goes here
    
    methods
        function obj = DOF_Stokes(filename,geometry) % Replace mesh for pdim
            nunkn_u = 2;
            nunkn_p = 1;
            obj.nunkn = [nunkn_u nunkn_p];
            [dirichlet_data,neumann_data,full_dirichlet_data] = Preprocess.getBC_fluids(filename,geometry);
            obj.getDOFconditions(geometry,dirichlet_data,neumann_data,full_dirichlet_data);
            obj.computeDOF(geometry);
        end
    end
end

