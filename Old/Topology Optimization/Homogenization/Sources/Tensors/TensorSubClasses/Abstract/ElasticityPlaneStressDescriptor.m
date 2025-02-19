classdef ElasticityPlaneStressDescriptor < ElasticityCaseDescriptor
    
    methods (Access = protected)
        
        function loadElasticityCaseVariable(obj)
            obj.elasticityCase = 'planeStress';
        end
     
    end
    
end

