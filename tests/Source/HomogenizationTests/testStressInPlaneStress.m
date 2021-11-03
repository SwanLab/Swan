classdef testStressInPlaneStress < handle
    
    properties (Access = public)
        tol = 1e-12;
    end
    
    properties (Access = private)
        stressFromTensorProduct
        stressFromVoigtProduct
        strain
        strainVoigtPS
        Ch
        ChVoigt
        ChVoigtPS
    end
    
    methods (Access = public)
        
        function obj = testStressInPlaneStress()
            obj.computeConstitutiveTensor()
            obj.computeStrain()
            obj.computeStressFromTensorProduct()
            obj.computeStressFromVoigtProduct()
        end
        
        function error = computeError(obj)
            sV  = obj.stressFromVoigtProduct;
            sT = obj.stressFromTensorProduct;
            error = norm(sV - sT);
        end

    end
    
    methods (Access = private)
        
        function computeConstitutiveTensor(obj)
            obj.Ch = Stiffness3DTensor;
            obj.Ch.createRandomTensor();
            obj.ChVoigt   = Tensor2VoigtConverter.convert(obj.Ch);
            obj.ChVoigtPS = PlaneStressTransformer.transform(obj.ChVoigt);
        end
        
        function computeStrain(obj)
            obj.strain = Strain3DTensor();
            obj.strain.makeItPlaneStressCompatible(obj.Ch);
            strainVoigt = Tensor2VoigtConverter.convert(obj.strain);
            obj.strainVoigtPS = PlaneStressTransformer.transform(strainVoigt);
        end
        
        function computeStressFromTensorProduct(obj)
            stress = ProductComputer.compute(obj.Ch,obj.strain);
            stressVoigt = Tensor2VoigtConverter.convert(stress);
            stressPS = PlaneStressTransformer.transform(stressVoigt);
            obj.stressFromTensorProduct = stressPS.getValue();
        end
        
        function computeStressFromVoigtProduct(obj)
            stressPS = ProductComputer.compute(obj.ChVoigtPS,obj.strainVoigtPS);
            obj.stressFromVoigtProduct = stressPS.getValue();
        end
        
    end

end