classdef SettingsCC < AbstractSettings
    
    properties (Access = protected, Abstract)
        defaultParamsName
    end
    
    properties (Access = public)
        shapeFuncSettings
        nShapeFuncs
        designVar
        homogenizedVarComputer
        targetParameters
        femData
        filterType
        physicalProblem
    end
    
    methods (Access = public)
        
        function obj = SettingsCC(varargin)
            if nargin == 1
                obj.loadParams(varargin{1});
            end
        end
        
        function list = getShapeFuncList(obj)
            list = obj.getShapeFunctionProperty('type');
        end
        
    end
    
    methods (Access = protected)
        
        function init(obj)
            obj.createShapeFunctionsSettings();
        end
        
        function list = getShapeFunctionProperty(obj,prop)
            nSF = obj.nShapeFuncs;
            list = cell(1,nSF);
            for iSF = 1:nSF
                list{iSF} = obj.shapeFuncSettings{iSF}.(prop);
            end
        end
        
    end
    
    methods (Access = private)
        
        function createShapeFunctionsSettings(obj)
            cParams = obj.shapeFuncSettings;
            nSF = length(cParams);
            sfS = cell(nSF,1);
            for iSF = 1:nSF
                if iscell(cParams)
                    s = cParams{iSF};
                else
                    s = cParams(iSF);
                end
                s.femSettings.fileName = obj.femData.fileName;
                s.femSettings.scale = obj.femData.scale;
                s.filterParams = obj.createFilterSettings();
                sfS{iSF} = SettingsShapeFunctional().create(s);
            end
            obj.shapeFuncSettings = sfS;
            obj.nShapeFuncs = nSF;
        end
        
        function s = createFilterSettings(obj)
            s.filterType  = obj.filterType;
            s.femSettings.fileName = obj.femData.fileName;
            s.femSettings.scale = obj.femData.scale;
            s = SettingsFilter(s);
        end
        
    end
    
end