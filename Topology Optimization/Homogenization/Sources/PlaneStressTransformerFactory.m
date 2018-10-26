classdef PlaneStressTransformerFactory < handle

    properties (Access = private)
        ClassName
        TensorCase
        PSTransformer
        Tensor
    end
    
    
    methods (Access = public)
        
        function PSTransformer = create(obj,Tensor)
            obj.init(Tensor)
            obj.createPlaneStressTransformer()
            PSTransformer = obj.getPlaneStressTransformer();
      end
    end
    
    methods (Access = private)
        
        function init(obj,Tensor)
            obj.Tensor = Tensor;
            obj.obtainClassName()
            obj.obtainTensorCase()
        end
        
        function obtainClassName(obj)
           obj.ClassName = class(obj.Tensor);
        end
        
        function obtainTensorCase(obj)
            if obj.isVoigt(obj.ClassName)
               [m,n] = size(obj.Tensor);
               if obj.isVoigtFourthOrderTensor(m,n)
                   obj.TensorCase = 'VoigtFourtOrder';
               elseif obj.isVoigtSecondOrderTensor(m,n)
                   obj.TensorCase = 'VoigtSecondOrder';                   
               else
                  error('Not admitted object to make it Plane Stress')
               end
            else
                obj.TensorCase = obj.ClassName;
            end  
        end
        
        function createPlaneStressTransformer(obj)
            switch obj.TensorCase
                case 'fourthOrderTensor'
                    obj.PSTransformer = PlaneStressTransformerForFourthOrderTensor(obj.Tensor);
                case 'secondOrderTensor'
                    obj.PSTransformer = PlaneStressTransformerForSecondOrderTensor(obj.Tensor);
                case {'VoigtFourtOrder'}
                    %obj.PSTransformer = PST4VoigtFourthOrderTensor(obj.Tensor);
                    obj.PSTransformer = PlaneStressVoigtTensorSymbolicallyTransformer(obj.Tensor);
                case {'VoigtSecondOrder'}
                    obj.PSTransformer = PST4VoigtSecondOrderTensor(obj.Tensor);
                otherwise
                    error('Not admitted object to make it Plane Stress')
            end            
        end
        
        function PS = getPlaneStressTransformer(obj)
            PS = obj.PSTransformer;
        end

    end

    methods (Static,Access = private)
        
        function ItIs = isVoigtFourthOrderTensor(m,n) 
           ItIs = m == 6 && n==6;
        end
        
        function ItIs = isVoigtSecondOrderTensor(m,n)
            ItIs = m == 1 && n == 6 || m ==6 && n == 1;
        end
        
        function isVoigt = isVoigt(ClassName)
            isDouble = strcmp(ClassName,'double');
            isSym    = strcmp(ClassName,'sym');
            isVoigt = isDouble || isSym;
        end
        
    end
    
    
    
    
    
end

