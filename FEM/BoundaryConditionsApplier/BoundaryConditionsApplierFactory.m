classdef BoundaryConditionsApplierFactory < handle
    
    methods (Access = public, Static)
        
        function obj = create(cParams)
            switch cParams.scale
                case 'MACRO'
                    switch cParams.type
                        case 'Neumann'
                            obj = NeumannConditionsApplier();
                        otherwise
                            obj = DirichletConditionsApplier(cParams);
                    end
                case 'MICRO'
                    obj = PeriodicBoundaryConditionApplier(cParams);
            end
        end
        
    end
    
end


