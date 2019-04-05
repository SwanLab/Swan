classdef SettingsShFunc_Chomog < SettingsShapeFunctional
  
    properties
        alpha
        beta
    end
    
     methods (Access = public)
        
        function obj = SettingsShFunc_Chomog(varargin)            
            switch nargin
                case 1
                    obj.loadParams(varargin{1});
            end
        end        
     end
end