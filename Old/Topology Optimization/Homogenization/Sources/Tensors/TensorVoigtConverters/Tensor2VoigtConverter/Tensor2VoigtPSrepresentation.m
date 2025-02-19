classdef Tensor2VoigtPSrepresentation < handle
    
    properties (Access = protected)
        indexTransformer
    end
    
    methods (Access = protected)
        
        function createIndexTransformer(obj)
            obj.indexTransformer = TensorVoigtIndexTransformerPS();
        end
        
    end
    
    
end