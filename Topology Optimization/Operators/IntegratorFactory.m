classdef IntegratorFactory
    methods (Access = public)
        function integrator = create(obj,mesh)
            type = obj.identifyMeshType(mesh);
            switch type
                case 'COMPOSITE'
                    integrator = Integrator_Composite(mesh);
                case 'INTERIOR'
                    mesh_composed = Mesh_Composed
                    integrator = Integrator_Interior(mesh_composed);
                case 'BOUNDARY'
                    integrator = Integrator_Boundary(mesh);
%                     integrator = Integrator_Unfitted(mesh);
                case 'GENERAL'
                    integrator = Integrator_General(mesh);
                otherwise
                    error('Invalid Mesh type. Currently, integrator only works with INTERIOR and BOUNDARY.')
            end
        end
    end
    
    methods (Static, Access = private)
        function type = identifyMeshType(mesh)
            if contains(class(mesh),'COMPOSITE','IgnoreCase',true)
                type = 'COMPOSITE';
            elseif contains(class(mesh),'BOUNDARY','IgnoreCase',true)
                type = 'BOUNDARY';
            elseif contains(class(mesh),'INTERIOR','IgnoreCase',true)
                type = 'INTERIOR';
            elseif contains(class(mesh),'MESH','IgnoreCase',true)
                type = 'GENERAL';
            end
        end
    end
end

