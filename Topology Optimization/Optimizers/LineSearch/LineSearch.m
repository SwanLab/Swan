classdef LineSearch < handle
    properties
        kappa
        kappa_min
        kfrac
    end
    
    methods (Static)
        function obj = create(settings,epsilon)
            switch settings.line_search
                case 'DIMENSIONALLY CONSISTENT'
                    switch settings.optimizer
                        case 'PROJECTED GRADIENT'
                            obj = LS_BackTracking_DimensionallyConsistent_PG(settings,epsilon);
                        case 'SLERP'
                            obj = LS_BackTracking_DimensionallyConsistent_SLERP;
                        case 'HAMILTON-JACOBI'
                            obj = LS_BackTracking_DimensionallyConsistent_HJ(settings);
                        otherwise
                            error('%s is NOT a valid unconstrained optimizer.',settings.optimizer);
                    end
                case 'DOUBLING LAST STEP'
                    switch settings.optimizer
                        case 'PROJECTED GRADIENT'
                            obj = LS_BackTracking_DoublingLastStep_PG(settings,epsilon);
                        case 'SLERP'
                            obj = LS_BackTracking_DoublingLastStep_SLERP;
                        case 'HAMILTON-JACOBI'
                            obj = LS_BackTracking_DoublingLastStep_HJ(settings);
                        otherwise
                            error('%s is NOT a valid unconstrained optimizer.',settings.optimizer);
                    end
                otherwise
                    error('%s is NOT a valid line-search algorithm.',settings.line_search);
            end
        end
    end
end

