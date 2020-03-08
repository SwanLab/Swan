classdef MeshPlotter_Null < MeshPlotter
    methods (Access = public, Static)
        function plot(~,~)
            warning('Cannot plot Volumetric meshes. Plot the corresponding Boundary instead.')
        end
    end
end