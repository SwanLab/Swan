classdef Mesh_Unfitted_3D_Interior < Mesh_Unfitted_3D & Mesh_Unfitted_Interior
    methods
        function obj = Mesh_Unfitted_3D_Interior(mesh_background,background_geom_interpolation)
            obj.storeBackgroundMesh(mesh_background,background_geom_interpolation);
            obj.max_subcells = 20;
            obj.nnodes_subcell = 4;
        end
        
        function V = computeVolume(obj)
            V = obj.computeMass;
        end
        
        function add2plot(~,~)
            warning('Solid bodies are not plotted.')
        end
        
        function computeDvoluCut(obj)
            x1 = obj.coord_iso_per_cell(:,1,1);   y1 = obj.coord_iso_per_cell(:,1,2);   z1= obj.coord_iso_per_cell(:,1,3);
            x2 = obj.coord_iso_per_cell(:,2,1);   y2 = obj.coord_iso_per_cell(:,2,2);   z2= obj.coord_iso_per_cell(:,2,3);
            x3 = obj.coord_iso_per_cell(:,3,1);   y3 = obj.coord_iso_per_cell(:,3,2);   z3= obj.coord_iso_per_cell(:,3,3);
            x4 = obj.coord_iso_per_cell(:,4,1);   y4 = obj.coord_iso_per_cell(:,4,2);   z4= obj.coord_iso_per_cell(:,4,3);
            
            J = x1.*y3.*z2-x1.*y2.*z3+x2.*y1.*z3-x2.*y3.*z1-x3.*y1.*z2+x3.*y2.*z1+x1.*y2.*z4-x1.*y4.*z2-x2.*y1.*z4+x2.*y4.*z1+...
                x4.*y1.*z2-x4.*y2.*z1-x1.*y3.*z4+x1.*y4.*z3+x3.*y1.*z4-x3.*y4.*z1-x4.*y1.*z3+x4.*y3.*z1+x2.*y3.*z4-x2.*y4.*z3...
                -x3.*y2.*z4+x3.*y4.*z2+x4.*y2.*z3-x4.*y3.*z2;
            obj.dvolu_cut = J/6;
        end
    end
end

