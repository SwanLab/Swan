classdef DOF_DiffReact < DOF
    
    methods (Access = public)
        
        function obj = DOF_DiffReact(mesh,interp) % Replace mesh for pdim
            obj.nunkn = 1;
            obj.dirichlet{1} = [];
            obj.dirichlet_values{1} = [];
            %% !! INITIALZE NEUMANN AS ALL NODES !!
            obj.neumann = [];
            obj.neumann_values  = [];
            obj.computeDOF(mesh,interp);
        end
        
    end
    
end

