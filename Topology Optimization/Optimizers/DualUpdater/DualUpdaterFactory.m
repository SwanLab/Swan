classdef DualUpdaterFactory < handle
    
   
    methods (Access = public, Static)
    
        function d = create(cParams)
            switch cParams.type
                case 'AlternatingPrimalDual'
                    d = DualUpdater_AugmentedLagrangian(cParams);
                case 'LagrangeMultiplierEstimation'
                    
                case 'NullSpace'
                    d = DualUpdater_NullSpace(cParams);
                case 'IPM'
                    d = DualUpdater_IPM(cParams);
            end
        end
    
    end
    
end