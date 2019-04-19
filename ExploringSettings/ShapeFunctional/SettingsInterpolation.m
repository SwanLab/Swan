classdef SettingsInterpolation < AbstractSettings
    
    
    properties (Access = protected)
        defaultParamsName = 'paramsMaterialInterpolation'
    end
   
    properties (Access = public)
        constitutiveProperties = struct
        typeOfMaterial
        interpolation 
        dim 
        type 
    end
    
    methods (Access = public)
        function obj = SettingsInterpolation(varargin)            
            if nargin == 1
                obj.loadParams(varargin{1})
            end
        end
    end
end