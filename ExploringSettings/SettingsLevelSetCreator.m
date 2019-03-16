classdef SettingsLevelSetCreator < DefaultSettings
 
    properties (Access = public)
        levelSetType = 'full'
        ndim = 2
        coord = Mesh_GiD('RVE_Square_Triangle_Fine.m').coord;
    end
    methods
    end
end