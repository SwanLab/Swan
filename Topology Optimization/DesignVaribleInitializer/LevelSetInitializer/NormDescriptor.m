classdef NormDescriptor < handle
    
    properties (Access = protected, Abstract)
        pos
        dist
    end
    
    methods (Access = protected, Abstract)
        computeDistance(obj)
    end
    
end

