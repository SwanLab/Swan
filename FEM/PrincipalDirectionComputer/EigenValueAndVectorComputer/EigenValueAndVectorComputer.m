classdef EigenValueAndVectorComputer < handle
    
    properties (Access = public)
        eigenVectorFunction
        eigenValueFunction
    end
    
    properties (Access = protected)
        ndim
    end
    
    methods (Access = public, Static)
        
        function obj = create(cParams)
            f = EigenValueAndVectorComputerFactory();
            obj = f.create(cParams);
        end
        
    end
    
    methods (Access = protected)
        
        function init(obj,cParams)
            obj.ndim = cParams.ndim;
        end
        
    end
    
end