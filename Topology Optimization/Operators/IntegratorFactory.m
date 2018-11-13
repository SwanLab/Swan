classdef IntegratorFactory
    methods (Access = public)
        function integrator = create(obj,mesh)
            type = obj.identifyMeshType(mesh);
            switch type
                case 'INTERIOR'
                    integrator = Integrator_Interior;
                case 'BOUNDARY'
                    integrator = Integrator_Boundary;
                otherwise
                    error('Invalid Mesh type. Currently, integrator only works with INTERIOR and BOUNDARY.')
            end
        end
    end
    
    methods (Static, Access = private)
        function type = identifyMeshType(mesh)
            if contains(class(mesh),'UNFITTED','IgnoreCase',true)
                if contains(class(mesh),'BOUNDARY','IgnoreCase',true)
                    type = 'BOUNDARY';
                else
                    type = 'INTERIOR';
                end
            else
                type = 'NOT UNFITTED';
            end
        end
    end
end

