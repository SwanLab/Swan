classdef SettingsHomogenizedVarComputerFromVademecum < ....
        AbstractSettings
    
    properties (Access = protected)
        defaultParamsName = 'paramsHomogenizedVarComputerFromVademecum'
    end
    
    properties (Access = public)
        fileName
    end
    
    methods (Access = public)
        function obj = SettingsHomogenizedVarComputerFromVademecum(varargin)
            if nargin == 1
                obj.loadParams(varargin{1})
            end
        end
    end
    
    
end