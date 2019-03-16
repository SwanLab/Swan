classdef Settings_TestLevelSetCreatorCustom < SettingsLevelSetCreator
 
    properties (Access = public)
    end
    
    methods (Access = public)
        
        function obj =Settings_TestLevelSetCreatorCustom()
            obj.coord = Mesh_GiD('test2d_quad.m').coord;
            obj.levelSetType = 'circle';
        end
        
    end
    
end