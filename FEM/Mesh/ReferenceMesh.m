classdef ReferenceMesh < handle
    
    properties (Access = public)
        coord
        connec
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = ReferenceMesh(cParams)
            obj.init(cParams)
        end
        
        function xGauss = computeXgauss(obj,xV)
            xGauss = xV;
        end

    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.coord = cParams.coord;
            obj.connec = cParams.connec;
        end
        
    end
    
end