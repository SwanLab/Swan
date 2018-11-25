classdef NormDescriptor < handle
    
    properties (Access = protected)
        dist
    end
    
    methods (Access = protected, Abstract)
        computeDistance(obj)
    end
    
end

