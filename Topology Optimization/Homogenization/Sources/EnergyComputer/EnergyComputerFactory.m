classdef EnergyComputerFactory < handle
    
    properties (Access = private)
        fourthOrder
        secondOrder
        energyComputer
    end
    
    methods (Access = public)
        
        function eComp = create(obj,fourthOrder,secondOrder)
            obj.init(fourthOrder,secondOrder)
            obj.createEnergyComputer()
            eComp = obj.getEnergyComputer();
        end
    end
    
    methods (Access = private)
       
        function init(obj,C,s)
            obj.fourthOrder = C;
            obj.secondOrder = s;
        end
        
        function createEnergyComputer(obj)
             C = obj.fourthOrder;
             s = obj.secondOrder;
             if obj.isVoigt()
                 e = EnergyComputerForVoigt(C,s);
             elseif obj.isTensor()
                  e = EnergyComputerForTensor(C,s);
             else
                 error('Not possible to compute Energy with these classes')
             end
             obj.energyComputer = e;
        end
        
        function e = getEnergyComputer(obj)
            e = obj.energyComputer;
        end
        
        function itIs = isVoigt(obj)
            isFourthVoigt = strcmp(obj.fourthOrder.getRepresentation(),'voigt');
            isSecondVoigt = strcmp(obj.secondOrder.getRepresentation(),'voigt');
            itIs = isFourthVoigt && isSecondVoigt;
        end
        
        function itIs = isTensor(obj)
            isFourthVoigt = strcmp(obj.fourthOrder.getRepresentation(),'tensor');
            isSecondVoigt = strcmp(obj.secondOrder.getRepresentation(),'tensor');
            itIs = isFourthVoigt && isSecondVoigt;
        end
        
    end
    
end

