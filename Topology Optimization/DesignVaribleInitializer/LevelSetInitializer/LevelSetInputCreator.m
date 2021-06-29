classdef LevelSetInputCreator < handle
    
    properties (Access = private)
        input
    end
    
    methods (Access = public)
        
        function obj = LevelSetInputCreator(settings,mesh)
            obj.input = input;
        end
        
        function i = getValue(obj)
            i = obj.input;            
        end
        
    end
    
    methods (Access = private, Static)
        

        
    end
    
    
    
end

