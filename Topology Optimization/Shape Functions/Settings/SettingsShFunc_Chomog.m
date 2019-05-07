classdef SettingsShFunc_Chomog < SettingsShapeFunctional
    
    properties (Access = public)
        alpha
        beta
    end
    
    methods (Access = public)
        
        function obj = SettingsShFunc_Chomog(varargin)
            obj.type = 'chomog_alphabeta';
            if nargin == 1
                obj.loadParams(varargin{1});
            end
        end
    end
end