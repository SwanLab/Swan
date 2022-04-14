classdef testComplianceTensorThrougtVoigtComparingEnergy < handle
    
    properties (Access = public)
         tol = 1e-10;
    end
    
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
           obj.computeEnergies();
       end

       function error = computeError(obj)
           es = obj.energyStiffTensProd;
           ec = obj.energyCompTensProd;
           error = norm(es - ec);
       end

    end
    
    methods (Access = private)
        
        function computeEnergies(obj)
            obj.generateStrain();
            obj.generateFourthOrderTensor();
            obj.computeStressTensor();
            obj.computeComplianceTensor();
            obj.computeEnergyByStifnessTensorProduct();
            obj.computeEnergyByComplianceTensorProduct();
        end
        
        function generateStrain(obj)
            obj.strain = Strain3DTensor;
            obj.strain.createRandomTensor();
        end
        
        function generateFourthOrderTensor(obj)
            obj.stiffTensor = Stiffness3DTensor;
            obj.stiffTensor.createRandomTensor();
        end
        
        function computeStressTensor(obj)
            e = obj.strain;
            c = obj.stiffTensor;
            s = ProductComputer.compute(c,e);
            obj.stress = s;
        end
        
        function computeComplianceTensor(obj)
            obj.compTensor = SymmetricFourthOrder3DTensor();
            c = obj.stiffTensor;
            invC = Inverter.invert(c);
            obj.compTensor = invC;
        end
        
        function computeEnergyByComplianceTensorProduct(obj)
            s = obj.stress;
            invc = obj.compTensor;
            energy = EnergyComputer.compute(invc,s);
            obj.energyStiffTensProd = energy;
        end
        
        function computeEnergyByStifnessTensorProduct(obj)
            e = obj.strain;
            c = obj.stiffTensor;
            energy = EnergyComputer.compute(c,e);
            obj.energyCompTensProd = energy;
        end
        
    end

end