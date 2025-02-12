classdef Voigt2TensorPSrepresentation < handle
    
    properties (Access = protected)
        indexTransformer
    end
    
    methods (Access = protected)
        
        function createIndexTransformer(obj)
            obj.indexTransformer = TensorVoigtIndexTransformerPS();
        end
        
    end
    
    
end