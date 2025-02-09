classdef PlaneStressTransformerFactory < handle

    properties (Access = private)
        PSTransformer
        tensor
    end
    
    
    methods (Access = public)
        
        function psTransformer = create(obj,Tensor)
            obj.init(Tensor)
            obj.createPlaneStressTransformer()
            psTransformer = obj.getPlaneStressTransformer();
      end
    end
    
    methods (Access = private)
        
        function init(obj,tensor)
            obj.tensor = tensor;
        end
        
        function createPlaneStressTransformer(obj)
            t = obj.tensor;
            ps = [];
            if obj.isVoigt()
                if obj.isStiffnessTensor()
                    if obj.isInvertible()
                        ps = PST4VoigtFourthOrderTensorNumerically(t.getValue());
                    else
                        ps = PST4VoigtFourthOrderTensorSymbolically(t.getValue());
                    end

                elseif obj.isSecondOrder()
                    ps = PST4VoigtSecondOrderTensor(t);
                end
                
            elseif obj.isTensor()
                if obj.isSecondOrder()
                    ps =  PlaneStressTransformerForSecondOrderTensor(t);
                end
                
            end
            
            if isempty(ps)
                error('Not admitted object to make it Plane Stress')
            else
                obj.PSTransformer = ps;                
            end

        end
        
        function PS = getPlaneStressTransformer(obj)
            PS = obj.PSTransformer;
        end
        
        function itIs = isVoigt(obj)
            itIs = strcmp(obj.tensor.getRepresentation(),'voigt');
        end
        
        function itIs = isTensor(obj)
            itIs = strcmp(obj.tensor.getRepresentation(),'tensor');
        end
        
        function ItIs = isSecondOrder(obj)
            ItIs = strcmp(obj.tensor.getOrder(),'second');
        end
        
         function itIs = isStiffnessTensor(obj)
            itIs = strcmp(obj.tensor.getFieldName,'stiffness');
        end
        
        function itIs = isInvertible(obj)
            t = obj.tensor.getValue();
            if obj.isSymbolic(t)
                itIs = isAlways(det(t) == 0);
            else
                itIs = min(abs(eig(t))) > 1e-13;
            end
        end
        
    end
    
    methods (Access = private, Static)
        function itIs = isSymbolic(a)
            itIs = isa(a,'sym');
        end
    end
    
    
    
end

