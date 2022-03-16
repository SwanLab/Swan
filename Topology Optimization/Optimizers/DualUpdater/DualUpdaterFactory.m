classdef DualUpdaterFactory < handle
    
   
    methods (Access = public, Static)
    
        function d = create(cParams)
            switch cParams.type
                case 'AugmentedLagrangian'
%                     d = DualUpdater_AugmentedLagrangian(cParams);
                    d = DualUpdater_NullSpace(cParams);
                case 'LagrangeMultiplierEstimation'

                case 'NullSpace'
                    d = DualUpdater_NullSpace(cParams);                   
            end
        end
    
    end
    
    methods (Access = private)
        
        
        
    end
   
    
end