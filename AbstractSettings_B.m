classdef AbstractSettings_B < handle
    
    properties (GetAccess = public, SetAccess = private)
        loadedFile
    end
    
    properties (Access = protected, Abstract)
        defaultParamsName
    end
    
    properties (Access = private)
        customParams
    end
    
    methods (Access = protected)
        
        function obj = AbstractSettings_B()
        end
        
    end
    
end