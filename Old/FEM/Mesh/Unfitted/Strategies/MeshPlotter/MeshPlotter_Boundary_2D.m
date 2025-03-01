classdef MeshPlotter_Boundary_2D < MeshPlotter
    methods (Access = public, Static)
        function plot(mesh,ax)
            hold on
            ncells = size(mesh.connec,1);
            for icell = 1:ncells
                plot(ax,mesh.coord(mesh.connec(icell,:),1),mesh.coord(mesh.connec(icell,:),2),'r.-','MarkerEdgeColor',[0.5 0 0],'MarkerSize',3);
            end
        end
    end
end