classdef ElasticityCaseDescriptor < handle
    
    properties (Access = protected)
        elasticityCase
    end
        
    methods (Access = public)
        
        function obj = ElasticityCaseDescriptor()
            obj.loadElasticityCaseVariable()
        end
        
        function e = getElasticityCase(obj)
            e = obj.elasticityCase;
        end

    end
    
    methods (Abstract,Access = protected)
        loadElasticityCaseVariable(obj)
    end
    
end
