classdef DOF_DiffReact < DOF
    %DOF_DiffReact Summary of this class goes here
    %   Detailed explanation goes here
    
    methods
        function obj = DOF_DiffReact(filename,geometry) % Replace mesh for pdim
            obj.nunkn = 1;
            [dirichlet_data,neumann_data,full_dirichlet_data] = Preprocess.getBC_mechanics(filename);
            obj.computeDOF(geometry,dirichlet_data,neumann_data,full_dirichlet_data);
        end
    end
end

