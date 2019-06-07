classdef BoundaryConditionsApplier < handle
    
    properties
    end
    
    methods (Access = public, Static)
        
        function bc = create(nfields,ndof,scale,type)
            switch type
                case 'Neumann'
                    bc = NeumannConditionsApplier();
                    %bc = DirichletConditionsApplier(nfields,ndof);                    
                otherwise
                    switch scale
                        case 'MACRO'
                            bc = DirichletConditionsApplier(nfields,ndof);
                        case 'MICRO'
                            bc = PeriodicBoundaryConditionApplier(nfields,ndof);
                    end
            end
            
        end
    end
    
end

