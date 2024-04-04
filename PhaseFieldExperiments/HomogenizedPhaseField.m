classdef HomogenizedPhaseField < handle
    
    properties (Access = private)
        fileName
        structuredMesh
        Ctensor
        microParams
    end
    
    methods (Access = public)
        
        function obj = HomogenizedPhaseField(cParams)
            obj.init(cParams)
            
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            
        end
        
    end
    
end