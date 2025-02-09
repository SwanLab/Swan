classdef VigdergauzParametersFactory < handle
    
    methods (Access = public, Static)

        function v = create(cParams)
            
           switch cParams.type               
               case 'VolumeAndRatio'
                   v = VigdergauzParametersFromVolumeAndPhi(cParams);                   
               case 'VolumeAndStrain'
                   v = VigdergauzParametersFromVolumeAndStrain(cParams);  
               case 'AxAndAy'
                   v = VigdergauzParametersFromAxAy(cParams);
           end
            
        end
        
    end
    
    
    
    
    
    
    
    
end