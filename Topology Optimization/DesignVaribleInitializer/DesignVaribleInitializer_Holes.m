classdef DesignVaribleInitializer_Holes < DesignVaribleInitializer
    properties
        N_holes
        R_holes
        phase_holes
    end
    
    methods
        function obj = DesignVaribleInitializer_Holes(settings,mesh,epsilon)
            obj@DesignVaribleInitializer(settings,mesh,epsilon);
            obj.load_holes_settings(settings);
        end
        
        function x = compute_initial_x(obj)
            phi = ones(size(obj.x));
            for i = 1:obj.mesh.ndim
                L(i) = max(obj.mesh.coord(:,i)) - min(obj.mesh.coord(:,i));
                phi = phi.*cos((obj.N_holes(i)+1)*(obj.mesh.coord(:,i)*pi)/L(i)+obj.phase_holes(i));
            end
            phi = phi + obj.R_holes-1;
            
            switch obj.optimizer
                case {'SLERP','HAMILTON-JACOBI'}
                    obj.x = phi;
                otherwise
                    initial_holes = ceil(max(phi,0))>0;
                    obj.x(initial_holes) = obj.hole_value;
            end
            
            bc = unique([obj.mesh.dirichlet(:,1); obj.mesh.pointload(:,1)]);
            if any(obj.x(bc)>0)
                warning('At least one BC is set on a hole')
            end
            x = obj.x;
        end
    end
    
    methods (Access = private)
        function load_holes_settings(obj,settings)
            obj.N_holes = settings.N_holes;
            obj.R_holes = settings.R_holes;
            obj.phase_holes = settings.phase_holes;
        end
    end
end

