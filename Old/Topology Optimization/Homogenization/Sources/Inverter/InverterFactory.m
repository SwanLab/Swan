classdef InverterFactory < handle
    
    properties (Access = private)
        tensor
        inverter
    end
    
    methods (Access = public)
        
        function inverter = create(obj,tensor)
            obj.init(tensor);
            obj.createInverter()
            inverter = obj.getInverter();
        end
    end

    methods (Access = private)
    
        function init(obj,tensor)
            obj.tensor = tensor;
        end

        function createInverter(obj)
            if obj.isVoigt()
                obj.inverter = VoigtTensorInverter(obj.tensor);
            elseif obj.isFourthOrderTensor() && ~obj.isSymmetric()
                obj.inverter = TensorInverter(obj.tensor); 
            elseif obj.isFourthOrderTensor() && obj.isSymmetric()
                obj.inverter = SymmetricTensorInverter(obj.tensor);                
            else
             error('Not admitted object to be Inverted')
            end
        end
        
        function inv = getInverter(obj)
            inv = obj.inverter;
        end
        
        function isVoigt = isVoigt(obj)
            isVoigt   = strcmp(obj.tensor.getRepresentation(),'voigt');
        end
        
        function itIs = isFourthOrderTensor(obj)
            itIs   = strcmp(obj.tensor.getOrder(),'fourth');
        end
        
        function isSymmetric = isSymmetric(obj)
            symmetrizer = FourthOrderSymmetrizer();
            isSymmetric = symmetrizer.isSymmetric(obj.tensor.getValue());            
        end

    end
    
end

