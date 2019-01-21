classdef MeshPlotter_Interior_3D < MeshPlotter_Abstract
    methods (Access = public, Static)
        function plot(~,~)
            warning('Cannot plot Volumetric meshes. Plot the corresponding Boundary instead.')
        end
    end
end