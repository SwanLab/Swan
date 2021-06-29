classdef SettingsHomogenizedVarComputerFromVademecum < ....
         SettingsHomogenizedVarComputer
    
    properties (Access = protected)
        defaultParamsName = 'paramsHomogenizedVarComputerFromVademecum.json'
    end
    
    properties (Access = public)
        fileName
        nelem
    end
    
    methods (Access = public)
        function obj = SettingsHomogenizedVarComputerFromVademecum(varargin)
            if nargin == 1
                obj.loadParams(varargin{1})
            end
        end
    end
    
    
end