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
            if ~isfield(varargin{1},'kFace')
                varargin{1}.kFace = 0;            
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
            s.nnode        = size(obj.connec,2);
            m = MeshTypeComputer(s);
            obj.type = m.compute();
        end
        
        function g = computeGeometryType(obj)
            sG.ndim           = size(obj.coord,2);
            sG.kFace          = obj.kFace;
            gC = GeometryTypeComputer(sG);            
            g = gC.compute();
        end
        
    end
    
end