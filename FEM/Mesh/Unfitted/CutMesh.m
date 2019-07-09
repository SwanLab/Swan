classdef CutMesh < handle
   
    properties (GetAccess = public, SetAccess = protected)
        coord
        connec
    end
        
    methods (Access = public)
        
        function obj = CutMesh(cParams)
            obj.coord  = cParams.coord;
            obj.connec = cParams.connec;
        end
        
    end         
    
end