classdef SettingsNumericalHomogenizer < AbstractSettings
    
    properties (Access = public)
        interpParams
        levelSetCreatorParams
        volumeShFuncParams
        elementDensityCreatorType
        outFileName
        testName
        print
        iter
        pdim
    end
    
    properties (Access = protected)
        defaultParamsName = 'paramsNumericalHomogenizer'
    end
    
    methods (Access = public)
        
        function obj = SettingsNumericalHomogenizer(varargin)
            if nargin == 1
                    obj.loadParams(varargin{1});
            end
            obj.setParams();
        end
        
    end
    
    methods (Access = private)
        
        function setParams(obj)
            obj.volumeShFuncParams.filename = obj.testName;
            obj.volumeShFuncParams.scale = 'MICRO';
        end
        
    end
    
end