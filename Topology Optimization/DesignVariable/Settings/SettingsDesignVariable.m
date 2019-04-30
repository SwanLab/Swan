classdef SettingsDesignVariable < AbstractSettings
    
    properties (Access = protected)
        defaultParamsName = 'paramsDesignVariable'
    end
    
    properties (Access = public)
        value
        mesh
        type
        initialCase
        levelSetCreatorSettings        
    end
    
     methods (Access = public)
        
        function obj = SettingsDesignVariable(varargin)
            if nargin == 1
                    obj.loadParams(varargin{1});
            end
        end 
        
     end
     
end