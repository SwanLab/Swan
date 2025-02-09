classdef testEnergyEquivalenceVoigtAndTensorNotation < handle
    
    properties (Access = public)
        tol = 1e-10;
    end

    properties (Access = protected)
        voigtEnergy
        tensorEnergy
        strain
        strainVoigt
        Ch
        ChVoigt
    end
    
    methods (Access = public)
        
        function obj = testEnergyEquivalenceVoigtAndTensorNotation()
            obj.init()
            obj.computeTensorEnergy()
            obj.computeVoigtEnergy();
        end
        
        function error = computeError(obj)
            error = abs(obj.voigtEnergy - obj.tensorEnergy);
        end
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.createStrains()
            obj.createConstitutiveTensors()
        end
        
        function createStrains(obj)
            obj.strain = Strain3DTensor();
            obj.strain.createRandomTensor();
            obj.strainVoigt = Tensor2VoigtConverter.convert(obj.strain);
        end
        
        function createConstitutiveTensors(obj)
            obj.generateFourthOrderTensor()
            obj.ChVoigt = Tensor2VoigtConverter.convert(obj.Ch);
        end
        
        function computeVoigtEnergy(obj)
            e = EnergyComputer.compute(obj.ChVoigt,obj.strainVoigt);
            obj.voigtEnergy = e;
        end
        
        function computeTensorEnergy(obj)
            e = EnergyComputer.compute(obj.Ch,obj.strain);
            obj.tensorEnergy = e;
        end
        
    end

    methods (Abstract)
       generateFourthOrderTensor(obj) 
    end

end