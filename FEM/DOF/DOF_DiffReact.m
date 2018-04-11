classdef DOF_DiffReact < DOF
    %DOF_DiffReact Summary of this class goes here
    %   Detailed explanation goes here    
    
    methods
        function obj = DOF_DiffReact(geometry) % Replace mesh for pdim
            obj.nunkn = 1;
            obj.dirichlet{1} = [];
            obj.dirichlet_values{1} = [];
            %% !! INITIALZE NEUMANN AS ALL NODES !!
            obj.neumann = [];
            obj.neumann_values  = [];
            obj.computeDOF(geometry);
        end
    end
end

