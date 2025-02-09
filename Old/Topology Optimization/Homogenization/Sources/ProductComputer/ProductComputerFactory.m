classdef ProductComputerFactory < handle
    
    properties (Access = private)
        secondOrder
        fourthOrder
        secondOrderOut
        productComputer
    end
    
    methods (Access = public)
        
        function pComp = create(obj,C,e)
            obj.init(C,e)
            obj.createProductComputer()
            pComp = obj.getProductComputer();
        end
    end
    
    methods (Access = private)
        
        function init(obj,C,e)
           obj.fourthOrder = C;
           obj.secondOrder = e;
        end
        
        function  createProductComputer(obj)
            C = obj.fourthOrder; 
            e = obj.secondOrder;
            
            if obj.isVoigt()
                p = ProductComputerForVoigt(C,e);
            elseif obj.isTensor()
                p = ProductComputerForTensor(C,e);
            else
                 error('Not possible to compute Product with these classes')
            end
            
            obj.productComputer = p;
            
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
        
        function p = getProductComputer(obj)
            p = obj.productComputer;
        end

        
    end
    
end

