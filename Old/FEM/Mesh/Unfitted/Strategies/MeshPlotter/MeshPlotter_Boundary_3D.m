classdef MeshPlotter_Boundary_3D < MeshPlotter
    methods (Access = public, Static)
        function plot(mesh,ax)
            hold on
            patch(ax,'vertices',mesh.coord,'faces',mesh.connec,...
                'edgecolor',[0.5 0 0], 'edgealpha',0.5,'edgelighting','flat',...
                'facecolor',[1 0 0],'facelighting','flat')
        end
    end
end