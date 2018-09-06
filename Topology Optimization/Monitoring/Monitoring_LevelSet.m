classdef Monitoring_LevelSet < Monitoring
    properties
        unfitted_mesh
        geometry
    end
    
    methods
        function obj = Monitoring_LevelSet(settings,mesh,monitoring_ON,plotting_ON)
            obj@Monitoring(settings,mesh,monitoring_ON,plotting_ON);
            obj.geometry= Geometry(obj.mesh,'LINEAR');
        end
    end
end

