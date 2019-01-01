classdef PostProcessDataBaseCreatorWithGaussData < PostProcessDataBaseCreator
    
    methods (Access = public)
               
        function obj = PostProcessDataBaseCreatorWithGaussData(d)
            obj.compute(d)
            obj.storeGaussDataInfo(d)
        end  
        
        function storeGaussDataInfo(obj,d)
            obj.data.ngaus = d.quad.ngaus;
            obj.data.posgp = d.quad.posgp';             
        end
        
    end

     
end