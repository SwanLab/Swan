classdef SettingsDesignVariable < AbstractSettings_B
    
    properties (Access = protected)
        defaultParamsName = 'paramsDesignVariable'
    end
    
    properties (Access = public)
        value
        mesh
        type
        initialCase
        levelSetCreatorSettings        
        scalarProductSettings        
    end
    
     methods (Access = public)
        
        function obj = SettingsDesignVariable()
        end 
        
     end
     
end