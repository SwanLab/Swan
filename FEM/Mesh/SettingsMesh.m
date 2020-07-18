classdef SettingsMesh < AbstractSettings
    
    properties (Access = protected)
        defaultParamsName = 'paramsMesh.json'
    end
    
    properties (GetAccess = public, SetAccess = public)
        coord
        connec
        type
        kFace
    end
    
    methods (Access = public)
        
        function obj = SettingsMesh(varargin)
            obj.loadParams(varargin{1});
            obj.computeType();
        end
        
    end
    
    methods (Access = private)
        
        function computeType(obj)
            s.ndim  = size(obj.coord,2);
            s.nnode = size(obj.connec,2);
            s.kFace = obj.kFace;
            m = MeshTypeComputer(s);
            obj.type = m.compute();
        end
        
    end
    
end