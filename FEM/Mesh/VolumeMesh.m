classdef VolumeMesh < Mesh
    
    properties (Access = public)
        geometryType = 'Volume';
    end
    
    properties (Access = private)
        cParams;
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = VolumeMesh(cParams)
            obj = obj@Mesh(cParams);
            obj.initVol(cParams)
            
        end
        
    end
    
    methods (Access = private)
        
        function initVol(obj,cParams)
        end
        
    end
    
end