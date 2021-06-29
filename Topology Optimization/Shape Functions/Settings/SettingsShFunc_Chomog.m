classdef SettingsShFunc_Chomog < SettingsShapeFunctional
    
    properties (Access = public)
        alpha
        beta
        ChTarget
    end
    
    methods (Access = public)
        
        function obj = SettingsShFunc_Chomog(varargin)
            obj.type = varargin{1}.type;
            if nargin == 1
                obj.loadParams(varargin{1});
            end
        end
    end
end