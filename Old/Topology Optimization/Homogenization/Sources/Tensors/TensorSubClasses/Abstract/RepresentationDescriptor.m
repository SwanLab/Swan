classdef RepresentationDescriptor < handle
    
    properties (Access = protected)
        representation
    end
        
    methods (Access = public)
        
        function obj = RepresentationDescriptor()
            obj.loadRepresentationVariable()
        end
        
        function r = getRepresentation(obj)
            r = obj.representation;
        end
    end
    
    methods (Abstract, Access = protected)
        loadRepresentationVariable(obj)
    end
    
end
