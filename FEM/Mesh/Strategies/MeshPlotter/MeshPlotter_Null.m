classdef MeshPlotter_Null < MeshPlotter_Abstract
    methods (Access = public, Static)
        function plot(~,~)
            warning('Cannot plot Volumetric meshes. Plot the corresponding Boundary instead.')
        end
    end
end