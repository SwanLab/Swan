classdef testComplianceTensorThrougtVoigtComparingEnergy < test
    
    properties (Access = private)
        strain
        stress
        stiffTensor
        compTensor
        energyStiffTensProd
        energyCompTensProd        
    end
    
    methods (Access = public)
        
        function obj = testComplianceTensorThrougtVoigtComparingEnergy()
            obj.computeEnergies
        end
    end
    
    methods (Access = private)
        
        function computeEnergies(obj)            
            obj.generateStrain()
            obj.generateFourthOrderTensor()
            obj.computeStressTensor()
            obj.computeComplianceTensor()
            obj.computeEnergyByStifnessTensorProduct()   
            obj.computeEnergyByComplianceTensorProduct()
        end
        
        function generateStrain(obj)
            obj.strain = StrainTensor();     
        end
        
        function generateFourthOrderTensor(obj)
            obj.stiffTensor = FourthOrderTensor;
            obj.stiffTensor.createRandomTensor();
        end
        
        function computeStressTensor(obj)
            producter = FourthWithSecondOrderProductComputer();
            e = obj.strain.getValue();
            c = obj.stiffTensor.getValue();
            s = producter.computeInTensor(c,e);
            obj.stress = StressTensor();
            obj.stress.setValue(s);
        end
        
        function computeComplianceTensor(obj)
            obj.compTensor = FourthOrderTensor();
            c = obj.stiffTensor;
            invC = Inverter.invert(c);
            obj.compTensor.setValue(invC);
        end
        
        function computeEnergyByComplianceTensorProduct(obj)
            s = obj.stress;
            invc = obj.compTensor;
            energy = obj.computeEnergy(s,invc);
            obj.energyStiffTensProd = energy;
        end
        
        function computeEnergyByStifnessTensorProduct(obj)
            e = obj.strain;
            c = obj.stiffTensor;
            energy = obj.computeEnergy(e,c);
            obj.energyCompTensProd = energy;
        end
        
    end
    
    methods (Access = private, Static)

        function en = computeEnergy(e,c)
            en = EnergyComputer.computeTensorEnergy(e,c);            
        end

    end
    
    methods (Access = protected)
       function hasPassed = hasPassed(obj)
           es = obj.energyStiffTensProd; 
           ec = obj.energyCompTensProd; 
           hasPassed = norm(es - ec) < 1e-12;
        end 
        
    end
    
end

