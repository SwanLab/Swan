classdef BoundaryConditionsApplierFactory < handle
    
    methods (Access = public, Static)
        
        function obj = create(cParams)
            switch cParams.scale
                case 'MACRO'
                    switch cParams.type
                        case 'Neumann'
                            % Obsolete
%                             obj = NeumannConditionsApplier();
                        otherwise
                            % Legacy
                            obj = DirichletConditionsApplier(cParams);
                    end
                case 'MICRO'
                    % Obsolete
%                     obj = PeriodicBoundaryConditionApplier(cParams);
            end
        end
        
    end
    
end
