classdef testAnisotropicPlaneStressbyEnergyEquivalence < test
    
    properties (Access = protected)
       
        
        VoigtEnergyInPlaneStress
        TensorEnergy
        
        stress
        strain
        Ch
        stress2PlaneStress
    end
    
    methods
        
        function obj = testAnisotropicPlaneStressbyEnergyEquivalence()
            obj.init()
            obj.computeVoigtEnergyInPlaneStress();
            obj.computeTensorEnergyInPlaneStress()
        end
        
        function init(obj)            
            obj.computeConstitutiveTensor()            
            obj.computeStrainTensor()
            obj.computeStressTensor()            
        end
        
        
        function computeStressTensor(obj)
            obj.stress = StressTensor();
            obj.stress.tensor      = obj.compute4and2orderTensorProduct(obj.Ch.tensor,obj.strain.tensor);
            obj.stress.makeItPlaneStress();
            obj.stress2PlaneStress = obj.compute4and2orderVoigtInPlaneStressProduct(obj.Ch.tensorVoigtInPlaneStress,obj.strain.tensorVoigtInPlaneStress);
        end
        
        function computeStrainTensor(obj)
            obj.strain = StrainTensor();
            obj.strain.createWithPlaneStressCompatibility(obj.Ch)
        end
        
        function computeConstitutiveTensor(obj)
            obj.Ch = fourthOrderTensor();
            obj.Ch.createRandomTensor();
            obj.Ch.computeTensorVoigt();
            obj.Ch.computeTensorVoigtInPlaneStress();
        end
        
        
        function computeVoigtEnergyInPlaneStress(obj)
            dim = size(obj.stress.tensorVoigtInPlaneStress,1);
            obj.VoigtEnergyInPlaneStress = 0;
            for i = 1:dim
                for j = 1:dim
                    StrainA = obj.strain.tensorVoigtInPlaneStress(i);
                    CH      = obj.Ch.tensorVoigtInPlaneStress(i,j);
                    StrainB = obj.strain.tensorVoigtInPlaneStress(j);
                    Energy = StrainA*CH*StrainB;
                    obj.VoigtEnergyInPlaneStress = obj.VoigtEnergyInPlaneStress + Energy;
                end
            end
        end
        
        function tensor = compute4and2orderTensorProduct(obj,fourthTensor,secondTensor)
            dim = size(secondTensor,1);
            tensor = zeros(size(secondTensor));
            for i = 1:dim
                for j = 1:dim
                    for k = 1:dim
                        for l = 1:dim
                            fourth = fourthTensor(i,j,k,l);
                            second = secondTensor(k,l);
                            tensor(i,j) = tensor(i,j) + fourth*second;
                        end
                    end
                end
            end
        end
        
        function tensor = compute4and2orderVoigtInPlaneStressProduct(obj,fourthTensor,secondTensor)
            dim = size(secondTensor,1);
            tensor = zeros(size(secondTensor));
            for i = 1:dim
                for j = 1:dim
                    fourth = fourthTensor(i,j);
                    second = secondTensor(j);
                    tensor(i) = tensor(i) + fourth*second;
                end
            end
        end
        
        function computeTensorEnergyInPlaneStress(obj)
            dim = size(obj.stress.tensor,1);
            obj.TensorEnergy = 0;
            for i = 1:dim
                for j = 1:dim
                    for k = 1:dim
                        for l = 1:dim
                            StrainA = obj.strain.tensor(i,j);
                            CH     = obj.Ch.tensor(i,j,k,l);
                            StrainB = obj.strain.tensor(k,l);
                            Energy = StrainA*CH*StrainB;
                            obj.TensorEnergy = obj.TensorEnergy + Energy;
                        end
                    end
                end
            end
        end
        
    end
    
    methods (Access = protected)
        
        function hasPassed = hasPassed(obj)
            isEqual1 = norm(double(obj.VoigtEnergyInPlaneStress) - obj.TensorEnergy) < 1e-6;
            isEqual2 = norm(double(obj.stress.tensorVoigtInPlaneStress) - obj.stress2PlaneStress(:)) < 1e-6;
            hasPassed = isEqual1 && isEqual2;
        end
        

    end
end

