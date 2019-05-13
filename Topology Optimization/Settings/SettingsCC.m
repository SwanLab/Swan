classdef SettingsCC < AbstractSettings
    
    properties (Access = protected, Abstract)
        defaultParamsName
    end
    
    properties (Access = public)
        settings
        shapeFuncSettings
        nShapeFuncs
        designVar
        homogenizedVarComputer
        targetParameters
        problemData
    end
    
    methods (Access = public)
        
        function obj = SettingsCC(varargin)
            if nargin == 1
                obj.loadParams(varargin{1});
            end
        end
        
    end
    
    methods (Access = protected)
        
        function init(obj)
            obj.createShapeFunctionsSettings();
        end
        
    end
    
    methods (Access = private)
        
        function createShapeFunctionsSettings(obj)
            s = obj.shapeFuncSettings;
            nSF = length(s);
            sfS = cell(nSF,1);
            for iSF = 1:nSF
                if iscell(s)
                    s{iSF}.filename = obj.problemData.problemFileName;
                    s{iSF}.scale = obj.problemData.scale;
                    s{iSF}.filterParams = obj.createFilterSettings();
                    sfS{iSF} = SettingsShapeFunctional().create(s{iSF},obj.settings);
                else
                    s(iSF).filename = obj.problemData.problemFileName;
                    s(iSF).scale = obj.problemData.scale;
                    s(iSF).filterParams = obj.createFilterSettings();
                    sfS{iSF} = SettingsShapeFunctional().create(s(iSF),obj.settings);
                end
            end
            obj.shapeFuncSettings = sfS;
            obj.nShapeFuncs = nSF;
        end
        
        function s = createFilterSettings(obj)
            s.filterType = obj.settings.filter;
            s = SettingsFilter(s);
        end
        
    end
    
end