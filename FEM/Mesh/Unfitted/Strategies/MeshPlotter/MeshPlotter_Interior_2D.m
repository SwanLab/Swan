classdef MeshPlotter_Interior_2D < MeshPlotter & PatchedMeshPlotter_Abstract
    methods (Access = protected, Static)
        function plotSurface(ax,coord,connec)
            hold on
            patch(ax,'vertices',coord,'faces',connec,...
                'edgecolor',[0.5 0 0], 'edgealpha',0.5,'edgelighting','flat',...
                'facecolor',[1 0 0],'facelighting','flat')
        end
    end
end