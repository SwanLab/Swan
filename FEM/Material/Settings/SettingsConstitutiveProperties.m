classdef SettingsConstitutiveProperties < AbstractSettings
    
    properties (Access = protected)
        defaultParamsName = 'paramsConstitutiveProperties.json'
    end
    
    properties (Access = public)
        rho_plus
        rho_minus
        E_plus
        E_minus
        nu_plus
        nu_minus
    end
    
    methods (Access = public)
        
        function obj = SettingsConstitutiveProperties(varargin)
            if nargin == 1
                obj.loadParams(varargin{1})
            end
        end
        
    end
       
end