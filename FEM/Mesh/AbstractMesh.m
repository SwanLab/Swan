classdef AbstractMesh < handle
    
    properties (GetAccess = public, SetAccess = protected)
        coord
        connec
        
        nelem
        geometryType
    end
    
end