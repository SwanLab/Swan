classdef DualUpdaterFactory < handle
    
   
    methods (Access = public, Static)
    
        function d = create(cParams)
            switch cParams.type
                case 'AugmentedLagrangian'
                    d = DualUpdater_AugmentedLagrangian(cParams);                    
                case 'LagrangeMultiplierEstimation'
                    
            end
        end
    
    end
    
    methods (Access = private)
        
        
        
    end
   
    
end