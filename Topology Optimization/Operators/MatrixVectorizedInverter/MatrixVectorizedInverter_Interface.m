classdef MatrixVectorizedInverter_Interface < handle
    
    methods (Access = public, Abstract)
        
        computeInverse(obj)
        computeDeterminant(obj)
        
    end
    
end