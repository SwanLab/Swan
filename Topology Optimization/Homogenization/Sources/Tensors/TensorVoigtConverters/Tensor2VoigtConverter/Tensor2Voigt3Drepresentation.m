classdef Tensor2Voigt3Drepresentation < handle
    
    properties (Access = protected)
        indexTransformer
    end
    
    methods (Access = protected)
        
        function createIndexTransformer(obj)
            obj.indexTransformer = TensorVoigtIndexTransformer3D();
        end
        
    end
    
    
end 
