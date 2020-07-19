classdef GeometryTypeComputer < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
        kFace
        ndim
    end
    
    methods (Access = public)
        
        function obj = GeometryTypeComputer(cParams)
            obj.init(cParams)            
        end
        
        function g = compute(obj)
            nGeom = obj.ndim + obj.kFace;
            switch nGeom
                case 1
                    g = 'Line';
                case 2
                    g = 'Surface';
                case 3
                    g = 'Volume';
            end            
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.ndim  = cParams.ndim;
            obj.kFace = cParams.kFace;
        end
        
    end
    
end