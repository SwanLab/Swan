classdef LineSearch < handle
    properties
        kappa
        kappa_min
        kfrac
    end
    
    methods (Static)
        function obj = create(cParams)
            switch cParams.type
                case 'DIMENSIONALLY CONSISTENT'
                    switch cParams.optimizerType
                        case 'PROJECTED GRADIENT'
                            obj = LS_BackTracking_DimensionallyConsistent_PG(cParams);
                        case 'SLERP'
                            obj = LS_BackTracking_DimensionallyConsistent_SLERP(cParams);
                        case 'HAMILTON-JACOBI'
                            obj = LS_BackTracking_DimensionallyConsistent_HJ(cParams);
                        otherwise
                            error('%s is NOT a valid unconstrained optimizer.',cParams.optimizer);
                    end
                case 'DOUBLING LAST STEP'
                    switch cParams.optimizerType
                        case 'PROJECTED GRADIENT'
                            obj = LS_BackTracking_DoublingLastStep_PG(cParams);
                        case 'SLERP'
                            obj = LS_BackTracking_DoublingLastStep_SLERP;
                        case 'HAMILTON-JACOBI'
                            obj = LS_BackTracking_DoublingLastStep_HJ(cParams);
                        otherwise
                            error('%s is NOT a valid unconstrained optimizer.',cParams.optimizer);
                    end
                otherwise  
                    error('%s is NOT a valid line-search algorithm.',cParams.type);
            end
        end
    end
end

