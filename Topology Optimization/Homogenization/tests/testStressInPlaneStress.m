classdef testStressInPlaneStress < test
    
    properties (Access = private)
        
        StressFromTensorProduct
        StressFromVoigtProduct
        
        strain
        Ch
        productComputer
    end
    
    methods (Access = public)
        
        function obj = testStressInPlaneStress()
            obj.productComputer = FourthWithSecondOrderProductComputer();
            obj.computeConstitutiveTensor()
            obj.computeStrain()
            obj.computeStressFromTensorProduct()
            obj.computeStressFromVoigtProduct()
        end
        
    end
    
    methods (Access = private)
        
        function computeConstitutiveTensor(obj)
            obj.Ch = FourthOrderTensor();
            obj.Ch.createRandomTensor();
            obj.Ch.computeTensorVoigt();
            obj.Ch.computeTensorVoigtInPlaneStress();
        end
        
        function computeStrain(obj)
            obj.strain = StrainTensor();
            obj.strain.createWithPlaneStressCompatibility(obj.Ch)
        end
        
        function computeStressFromTensorProduct(obj)
            Ctensor = obj.Ch.tensor;
            Strain  = obj.strain.tensor;
            obj.StressFromTensorProduct = StressTensor();
            Stress = obj.productComputer.computeInTensor(Ctensor,Strain);
            obj.StressFromTensorProduct.tensor = Stress;
            obj.StressFromTensorProduct.makeItPlaneStress();
        end
        
        function computeStressFromVoigtProduct(obj)
            Ctensor = obj.Ch.tensorVoigtInPlaneStress;
            Strain = obj.strain.tensorVoigtInPlaneStress;
            obj.StressFromVoigtProduct = StressTensor();
            Stress = obj.productComputer.computeInVoigt(Ctensor,Strain);
            obj.StressFromVoigtProduct.tensorVoigtInPlaneStress = Stress;
        end
        
    end
    
    methods (Access = protected)
        
        function hasPassed = hasPassed(obj)
            StFromVoigt  = obj.StressFromVoigtProduct.tensorVoigtInPlaneStress;
            StFromTensor = double(obj.StressFromTensorProduct.tensorVoigtInPlaneStress);
            hasPassed = norm(StFromVoigt - StFromTensor) < 1e-12;
        end
        
        
    end
end

