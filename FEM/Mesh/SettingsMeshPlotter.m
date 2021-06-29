classdef SettingsMeshPlotter < AbstractSettings
    
    properties (Access = protected)
        defaultParamsName = 'paramsMeshPlotter.json'
    end
    
    properties (Access = public)
        mesh
        isBackground
    end
    
    methods (Access = public)
        
        function obj = SettingsMeshPlotter(varargin)
            switch nargin
                case 1
                    obj.loadParams(varargin{1});
            end
        end
        
    end
        
end