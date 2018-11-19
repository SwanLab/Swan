classdef Mesh_Unfitted_2D_Interior < Mesh_Unfitted_2D & Mesh_Unfitted_Interior
    methods (Access = public)
        function obj = Mesh_Unfitted_2D_Interior(mesh_background,background_geom_interpolation)
            obj.storeBackgroundMesh(mesh_background,background_geom_interpolation);
            obj.max_subcells = 6;
            obj.nnodes_subcell = 3;
        end
        
        function add2plot(obj,h)
            hold on;
            patch(h,'vertices',obj.coord,'faces',obj.connec,...
                'edgecolor',[0.5 0 0], 'edgealpha',0.5,'edgelighting','flat',...
                'facecolor',[1 0 0],'facelighting','flat')
        end
        
        function S = computeSurface(obj)
            S = obj.computeMass;
        end
        
        function computeDvoluCut(obj)
            x1 = obj.coord_iso_per_cell(:,1,1);   y1 = obj.coord_iso_per_cell(:,1,2);
            x2 = obj.coord_iso_per_cell(:,2,1);   y2 = obj.coord_iso_per_cell(:,2,2);
            x3 = obj.coord_iso_per_cell(:,3,1);   y3 = obj.coord_iso_per_cell(:,3,2);
            obj.dvolu_cut = 0.5*abs((x2-x1).*(y3-y1)-(x3-x1).*(y2-y1));
        end
    end
end

