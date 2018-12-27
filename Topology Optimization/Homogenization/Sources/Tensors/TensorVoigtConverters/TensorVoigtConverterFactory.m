classdef TensorVoigtConverterFactory < handle
    
    
    properties (Access = protected)
        input
    end
    
    
    methods (Access = protected)
        
        function itIs = isStrainTensor(obj)
            itIs = strcmp(obj.input.getFieldName(),'strain');
        end
        
        function itIs = isStressTensor(obj)
            itIs = strcmp(obj.input.getFieldName(),'stress');
        end
        
        function itIs = isComplianceTensor(obj)
            itIs = strcmp(obj.input.getFieldName,'compliance');
        end

        function itIs = isStiffnessTensor(obj)
            itIs = strcmp(obj.input.getFieldName,'stiffness');
        end
        
        function itIs = is3D(obj)
            itIs = strcmp(obj.input.getElasticityCase,'3D');
        end
        
        function itIs = isPlaneStress(obj)
            itIs = strcmp(obj.input.getElasticityCase,'planeStress');
        end
        
        
    end
    
    methods (Access = private, Static)
        
        function showError()
            error('Not admitted object to transform')
        end
    end

    
end