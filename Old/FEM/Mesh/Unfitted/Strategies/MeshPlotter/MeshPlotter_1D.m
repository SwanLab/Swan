classdef MeshPlotter_1D < MeshPlotter & PatchedMeshPlotter_Abstract
    methods (Access = protected, Static)
        function plotSurface(ax,coord,connec)
            hold on
            ncells = size(connec,1);
            for icell = 1:ncells
                plot(ax,coord(connec(icell,:),1),coord(connec(icell,:),2),'r.-','MarkerEdgeColor',[0.5 0 0],'MarkerSize',3);
            end
        end
    end
end