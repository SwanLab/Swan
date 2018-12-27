classdef BoundaryConditionsApplier < handle
    
    properties
    end
    
    methods (Access = public, Static)
        
        function bc = create(nfields,ndof,scale)
            
            switch scale
                case 'MACRO'
                   bc = DirichletConditionsApplier(nfields,ndof);
                case 'MICRO'
                   bc = PeriodicBoundaryConditionApplier(nfields,ndof);                    
            end                                            
            
        end
    end
    
end

