classdef testEnergyEquivalenceVoigtAndTensorNotationForIsoTensor < test
    
    properties (Access = protected)
        VoigtEnergy
        TensorEnergy
        
        stress
        strain
        Ch
    end
    
    methods
        
        function obj = testEnergyEquivalenceVoigtAndTensorNotationForIsoTensor()
            obj.init()
            obj.computeVoigtEnergy();
            obj.computeTensorEnergy()
        end
        
        function init(obj)
            obj.strain = StrainTensor();
            obj.Ch = IsotropicConstitutiveTensor3D(1,1/3);
            obj.stress = StressTensor();
            obj.stress.tensor      = obj.compute4and2orderTensorProduct(obj.Ch.tensor,obj.strain.tensor);
            obj.stress.tensorVoigt = obj.compute4and2orderVoigtProduct(obj.Ch.tensorVoigt,obj.strain.tensorVoigt);
        end
        
        
        function computeVoigtEnergy(obj)
            dim = size(obj.stress.tensorVoigt,1);
            obj.VoigtEnergy = 0;
            for i = 1:dim
                for j = 1:dim
                    StrainA = obj.strain.tensorVoigt(i);
                    CH      = obj.Ch.tensorVoigt(i,j);
                    StrainB = obj.strain.tensorVoigt(j);
                    Energy = StrainA*CH*StrainB;
                    obj.VoigtEnergy = obj.VoigtEnergy + Energy;
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
        
        function tensor = compute4and2orderVoigtProduct(obj,fourthTensor,secondTensor)
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
        
        function computeTensorEnergy(obj)
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
            hasPassed = norm(double(obj.VoigtEnergy) - obj.TensorEnergy) < 1e-6;
        end


    end
end

