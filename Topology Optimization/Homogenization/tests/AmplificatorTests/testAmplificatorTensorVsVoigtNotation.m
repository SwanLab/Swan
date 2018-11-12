classdef testAmplificatorTensorVsVoigtNotation < test
    
    
    properties (Access = private)
        ampTensorVA
        ampTensorTA
        voigtAmpTensorTA
    end
       
    
    methods (Access = public)
        
        function obj = testAmplificatorTensorVsVoigtNotation() 
            obj.createConstitutiveTensors()
            obj.createCorrectorTensors()
            obj.computeVoigtAmplificatorByVoigtAlgebra()
            obj.computeVoigtAmplificatorByTensorAlgebra()
            obj.transformAmplificatorByTensorAlgebraInVoigt()
        end
      
        
    end
    
    methods (Access = private)
                
        
        function createConstitutiveTensors(obj) 
            obj.C = 
        end
        
        function createCorrectorTensors(obj)
        end
        
        function computeVoigtAmplificatorByVoigtAlgebra(obj)
            obj.ampTensorVA = SymmetricFourthOrder3DVoigtTensor();
            obj.ampTensorVA.createRandomTensor();
        end
        
        function computeVoigtAmplificatorByTensorAlgebra(obj)
            obj.ampTensorTA = SymmetricFourthOrder3DTensor();
            obj.ampTensorTA.createRandomTensor()
        end
        
        function transformAmplificatorByTensorAlgebraInVoigt(obj)
            ta = obj.ampTensorTA;
            obj.voigtAmpTensorTA = Tensor2VoigtConverter.convert(ta);
        end
        
    end
    
    
    methods (Access = protected)
        
        function hasPassed = hasPassed(obj)
            va = obj.ampTensorVA.getValue();
            ta = obj.voigtAmpTensorTA.getValue();
            hasPassed = norm(ta(:) - va(:))/norm(ta(:)) < 1e-12;
        end
        
    end
    

end

