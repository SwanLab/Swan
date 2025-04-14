classdef SettingsMesh < AbstractSettings
    
    properties (Access = protected)
        defaultParamsName% = 'paramsMesh.json'
    end
    
    properties (GetAccess = public, SetAccess = public)
        coord
        connec
        type
        kFace
        geometryType
        interpType
    end
    
    methods (Access = public)
        
        function obj = SettingsMesh(varargin)
            if ~max(isfield(varargin{1},'kFace'))
                varargin{1}.kFace = 0;
            end
            if (isfield(varargin{1},'interpType'))
                obj.interpType = varargin{1}.interpType;
            else                
                obj.interpType = 'LINEAR';
            end            
            obj.coord  = varargin{1}.coord;
            obj.connec = varargin{1}.connec;
            obj.kFace  = varargin{1}.kFace;
            %obj.loadParams(varargin{1});
            obj.computeType();
        end
        
    end
    
    methods (Access = private)
        
        function computeType(obj)
            s.geometryType = obj.computeGeometryType();
            s.nnodeElem    = size(obj.connec,2);
            m = MeshTypeComputer(s);
            obj.type = m.compute();
        end
        
        function g = computeGeometryType(obj)
            ndim  = size(obj.coord,2);
            nGeom = ndim + obj.kFace;
            switch nGeom
                case 1
                    g = 'Line';
                case 2
                    g = 'Surface';
                case 3
                    g = 'Volume';
            end
            obj.geometryType = g;
        end
        
    end
    
end