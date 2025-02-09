classdef Elasticity3dDescriptor < ElasticityCaseDescriptor
    
    methods (Access = protected)
        
        function loadElasticityCaseVariable(obj)
            obj.elasticityCase = '3D';
        end
     
    end
    
end

