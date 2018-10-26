classdef testAnisotropicPlaneStressbyEnergyEquivalence < test
    
    properties (Access = private)
        VoigtEnergyInPlaneStress
        TensorEnergy
        
        strain
        Ch
    end
    
    methods (Access = public)
        
        function obj = testAnisotropicPlaneStressbyEnergyEquivalence()
            obj.init()
            obj.computeVoigtEnergyInPlaneStress();
            obj.computeTensorEnergyInPlaneStress()
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)            
            obj.computeConstitutiveTensor()            
            obj.computeStrainTensor()
        end
       
        function computeStrainTensor(obj)
            obj.strain = StrainTensor();
            obj.strain.createWithPlaneStressCompatibility(obj.Ch)
        end
        
        function computeConstitutiveTensor(obj)
            obj.Ch = FourthOrderTensor();
            obj.Ch.createRandomTensor();
            obj.Ch.computeTensorVoigt();
            obj.Ch.computeTensorVoigtInPlaneStress();
        end
        
        function computeVoigtEnergyInPlaneStress(obj)
            EnComputer =  EnergyComputer(obj.strain,obj.Ch);
            obj.VoigtEnergyInPlaneStress = EnComputer.EnergyVoigtPlaneStress;
        end
                
        function computeTensorEnergyInPlaneStress(obj)
            En  = EnergyComputer.computeTensorEnergy(obj.strain,obj.Ch);
            obj.TensorEnergy = En;
        end
        
    end
    
    methods (Access = protected)
        
        function hasPassed = hasPassed(obj)
            EnVoigt = double(obj.VoigtEnergyInPlaneStress);
            EnTens  = obj.TensorEnergy;
            hasPassed = norm(EnVoigt - EnTens) < 1e-12;
        end
        

    end
end

