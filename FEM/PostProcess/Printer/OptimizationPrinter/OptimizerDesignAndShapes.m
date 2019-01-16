classdef OptimizerDesignAndShapes < OptimizerPrinterWithShapes
                                 
    methods (Access = protected) 
        
        function  createTopOptFields(obj,x,cost,constraint)
            obj.createTopOptFields@OptimizerPrinter(x,cost,constraint);            
            obj.createTopOptFields@OptimizerPrinterWithShapes(x,cost,constraint);            
        end 
        
        
    end
                                                           
end