classdef PostProcessDataBaseCreatorFactory < handle
    
    
    methods (Access = public, Static)

        function ps = create(hasGaussData,dI)
            
            if hasGaussData
                ps = PostProcessDataBaseCreatorWithGaussData(dI);
            else
                ps = PostProcessDataBaseCreatorWithNoGaussData(dI);
            end            
        end
      
    end 
    
    
end