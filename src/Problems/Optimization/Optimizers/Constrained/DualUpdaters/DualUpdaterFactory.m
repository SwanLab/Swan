classdef DualUpdaterFactory < handle
    
   
    methods (Access = public, Static)
    
        function d = create(cParams)
            switch cParams.type
                case 'Augmented Lagrangian'
                    d = DualUpdaterAugmentedLagrangian(cParams);
                case 'LagrangeMultiplierEstimation'
                    
                case 'NullSpace'
                    d = DualUpdaterNullSpace(cParams);
                case 'IPM'
                    d = DualUpdaterIPM(cParams);
            end
        end
    
    end
    
end