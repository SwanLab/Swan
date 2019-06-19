classdef AbstractMesh < handle
    
    properties (Access = public)
        unfittedType        
    end
    
    
    properties (GetAccess = public, SetAccess = protected)
        coord
        connec
        
        nelem
        geometryType
    end
    
end