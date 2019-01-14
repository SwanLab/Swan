classdef MeshPlotter_Interior_3D < MeshPlotter_Abstract
    methods (Access = public, Static)
        function plot(~,~)
            error('Cannot plot Volumetric meshes. Plot the corresponding Boundary instead.')
        end
    end
end