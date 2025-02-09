classdef VoigtRepresentation < RepresentationDescriptor
    
    properties (Access = protected)
        dimVoigt
    end 
    
    methods (Access = public)
        
        function d = getVoigtDimension(obj)
          d = obj.dimVoigt;
        end
    end
    
    methods (Access = protected)
        
        function obj = VoigtRepresentation()
            obj.loadVoigtDimensionValue()
        end
        
        function loadRepresentationVariable(obj)
            obj.representation = 'voigt';
        end
    end
    
    methods (Access = private)
        
        function loadVoigtDimensionValue(obj)
            tensSize = obj.getTensorSize();
            obj.dimVoigt = tensSize(1);
        end
        
    end
    
end

