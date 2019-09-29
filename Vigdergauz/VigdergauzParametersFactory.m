classdef VigdergauzParametersFactory < handle
    
    methods (Access = public, Static)

        function v = create(cParams)
            
           switch cParams.type               
               case 'VolumeAndRatio'
                   v = VigdergauzParametersFromThetaAndPhi(cParams);                   
               case 'VolumeAndStrain'
               
           end
            
        end
        
    end
    
    
    
    
    
    
    
    
end