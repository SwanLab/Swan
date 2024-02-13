classdef LineMesh < NewMesh
    
    properties (Access = public)
        geometryType = 'Line';
        
        coord, connec
        kFace
    end
    
    properties (Access = private)
        type
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = LineMesh(cParams)
            obj.init(cParams)
            obj.computeType();
            obj.createInterpolation();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.coord  = cParams.coord;
            obj.connec = cParams.connec;
            obj.kFace  = 0;
        end

        function computeType(obj)
%             s.geometryType = obj.geometryType;
%             s.nnodeElem    = obj.nnodeElem;
%             t = MeshTypeComputer(s);
%             obj.type = t.compute();
            obj.type = 'LINE';
        end

        function createInterpolation(obj)
            obj.interpolation = Interpolation.create(obj.type,'LINEAR');
        end
        
    end
    
end