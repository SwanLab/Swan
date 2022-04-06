classdef testAnisotropicPlaneStressbyEnergyEquivalence < handle
    
    properties (Access = public)
        tol = 1e-10;
    end
    
    properties (Access = private)
        energyFromPS
        energyFromTensor
        
        strain
        strainVoigt
        strainVoigtPS
        
        Ch
        ChVoigt
        ChVoigtPS
    end
    
    methods (Access = public)
        
        function obj = testAnisotropicPlaneStressbyEnergyEquivalence()
            obj.init()
            obj.computeVoigtEnergyInPlaneStress();
            obj.computeTensorEnergyInPlaneStress()
        end
        
        function error = computeError(obj)
            enPS = double(obj.energyFromPS);
            enTens  = obj.energyFromTensor;
            error =  norm(enPS - enTens);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.computeConstitutiveTensor()
            obj.computeStrainTensor()
        end
        
        
        function computeConstitutiveTensor(obj)
            obj.createChTensor()
            obj.createChVoigtTensor()
            obj.createChPlaneStressVoigtTensor()
        end
        
        function createChTensor(obj)
            obj.Ch = Stiffness3DTensor;
            obj.Ch.createRandomTensor();
        end
        
        function createChVoigtTensor(obj)
            c = Tensor2VoigtConverter.convert(obj.Ch);
            obj.ChVoigt = c;
        end
        
        function createChPlaneStressVoigtTensor(obj)
            c = PlaneStressTransformer.transform(obj.ChVoigt);
            obj.ChVoigtPS = c;
        end
        
        function computeStrainTensor(obj)
            obj.strain        = Strain3DTensor();
            obj.strain.makeItPlaneStressCompatible(obj.Ch);
            obj.strainVoigt   = Tensor2VoigtConverter.convert(obj.strain);
            obj.strainVoigtPS = PlaneStressTransformer.transform(obj.strainVoigt);
        end
        
        function computeVoigtEnergyInPlaneStress(obj)
            obj.energyFromPS = EnergyComputer.compute(obj.ChVoigtPS,obj.strainVoigtPS);
        end
        
        function computeTensorEnergyInPlaneStress(obj)
            obj.energyFromTensor = EnergyComputer.compute(obj.Ch,obj.strain);
        end
        
    end

end