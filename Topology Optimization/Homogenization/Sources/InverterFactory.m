classdef InverterFactory < handle
    
    properties (Access = private)
        tensor
        className
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
            obj.tensor = tensor.getValue();
            obj.obtainClassName()
        end
        
        function obtainClassName(obj)
            obj.className = class(obj.tensor);
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
            isDouble  = strcmp(obj.className,'double');
            isSym     = strcmp(obj.className,'sym');
            isTwoDim  = ismatrix(obj.tensor);
            isVoigt   = (isDouble || isSym ) && isTwoDim;
        end
        
        function itIs = isFourthOrderTensor(obj)
            %first  = isa(obj.ClassName,'fourthOrderTensor');
            %second = isa(obj.ClassName,'IsotropicConstitutiveTensor3D');
            isDouble   = strcmp(obj.className,'double');
            isSym      = strcmp(obj.className,'sym');
            isFourDim  = ndims(obj.tensor) == 4;
            itIs  = (isDouble || isSym ) && isFourDim;
        end
        
        function isSymmetric = isSymmetric(obj)
            symmetrizer = FourthOrderSymmetrizer();
            isSymmetric = symmetrizer.isSymmetric(obj.tensor);            
        end

    end
    
end

