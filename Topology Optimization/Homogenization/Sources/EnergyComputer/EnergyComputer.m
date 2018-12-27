classdef EnergyComputer < handle
    
    properties (Access = protected)
        energy
        secondOrder
        fourthOrder
    end
    
    methods (Access = public, Static)
        
        function e = compute(C,s)
            factory        = EnergyComputerFactory();
            energyComputer = factory.create(C,s);
            e = energyComputer.getEnergy();
        end
        
    end
    
    methods (Access = private)
        
        function e = getEnergy(obj)
            e = obj.energy;
        end
        
    end
    
    methods (Access = protected)
        
        function init(obj,C,s)
            obj.fourthOrder = C;
            obj.secondOrder = s;
        end
        
    end
    
    methods (Access = protected, Abstract)
        computeEnergy(obj)
    end
end

