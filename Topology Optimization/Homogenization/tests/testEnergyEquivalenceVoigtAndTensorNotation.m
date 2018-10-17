classdef testEnergyEquivalenceVoigtAndTensorNotation < test
    
    properties (Access = protected)
        VoigtEnergy
        TensorEnergy
        
        stress
        strain
        Ch
    end
    
    methods (Access = protected)
        
        function obj = testEnergyEquivalenceVoigtAndTensorNotation()
            obj.init()
            obj.computeVoigtEnergy();
            obj.computeTensorEnergy()
        end
        
        function hasPassed = hasPassed(obj)
            hasPassed = norm( double(obj.VoigtEnergy) - double(obj.TensorEnergy)) < 1e-6;
        end
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.strain = StrainTensor();
            obj.generateFourthOrderTensor();
        end
        
        function computeVoigtEnergy(obj)
            EnComputer =  EnergyComputer(obj.strain,obj.Ch);
            obj.VoigtEnergy = EnComputer.EnergyVoigt;
        end
        
        function computeTensorEnergy(obj)

            obj.TensorEnergy = EnergyComputer.computeTensorEnergy(obj.strain, ...
                                                                  obj.Ch);
        end
        
    end
    
    methods (Abstract)
       generateFourthOrderTensor(obj) 
    end
    

end

