classdef MeshPlotter_Boundary_2D < MeshPlotter_Abstract
    methods (Access = public, Static)
        function plot(mesh,ax)
            hold on
            ncells = size(mesh.connec,1);
            for icell = 1:ncells
                plot(ax,mesh.coord(mesh.connec(icell,:),1),mesh.coord(mesh.connec(icell,:),2),'k-');
            end
        end
    end
end