classdef TensorRepresentation < RepresentationDescriptor
    
    properties (Access = protected)
        dim 
    end
    
    methods (Access = public)
        
        function d = getDimension(obj)
            d = obj.dim;
        end
    end
    
    methods (Access = protected)
        
        function obj = TensorRepresentation()
            obj.loadDimensionValue()
        end
        
        function loadRepresentationVariable(obj)
            obj.representation = 'tensor';
        end
    end
    
    methods (Access = private)
        
        function loadDimensionValue(obj)
            tensSize = obj.getTensorSize();
            obj.dim = tensSize(1);
        end
        
    end
end

