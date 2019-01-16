classdef OptimizerElemntalDensityCase < OptimizerPrinterWithGaussData ...
                                  & OptimizerElementalDensity
                                 
    methods (Access = protected) 
        
        function  createTopOptFields(obj,x,cost,constraint)
            obj.createTopOptFields@OptimizerElementalDensity(x,cost,constraint);            
        end 
        
        
    end
                                                           
end