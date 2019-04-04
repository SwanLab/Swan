classdef SettingsFilterFactory < AbstractSettings
    
    properties (Access = protected)
        defaultParamsName = 'paramsFilterFactory'
    end
    
    properties (Access = public)
        type
        designVar
        domainType
    end
    
    methods (Access = public)
        
        function obj = SettingsFilterFactory(varargin)
            switch nargin
                case 1
                    obj.loadParams(varargin{1});
                case 2
                    obj.type = varargin{1};
                    obj.designVar = obj.getDesignVariableType(varargin{2});
            end
        end
        
    end
    
    methods (Access = private, Static)
        
        function designVar = getDesignVariableType(optimizer)
            switch optimizer
                case {'MMA','PROJECTED GRADIENT','IPOPT'}
                    designVar = 'DENSITY';
                case {'SLERP','HAMILTON-JACOBI','PROJECTED SLERP'}
                    designVar = 'LEVELSET';
            end
        end
        
    end
    
end