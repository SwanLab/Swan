classdef SettingsMeshPlotter < AbstractSettings
    
    properties (Access = protected)
        defaultParamsName = 'paramsMeshPlotter.json'
    end
    
    properties (Access = public)
        mesh
        isBackground
        faceColor
        faceAlpha
        edgeAlpha
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