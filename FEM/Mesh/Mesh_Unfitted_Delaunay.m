classdef Mesh_Unfitted_Delaunay < Mesh_Unfitted
    properties
    end
    methods
        function obj = computeCutMesh(obj)
            [nodes_n_cutpoints_iso,active_nodes] = obj.findCutPoints_Iso;
            nodes_n_cutpoints_global = obj.findCutPoints_Global;
            
            obj.computeDelaunay_allocateMemory;
            obj.i_debug=1;
            
            k0 = 0; m0 = 0; c0 = 0;
            k1 = 0; m1 = 0; c1 = 0;
            for icut = 1:length(obj.cut_cells)
                icell = obj.cut_cells(icut);
                subcell_cutPoints_iso = nodes_n_cutpoints_iso(active_nodes(:,:,icut),:,icut);
                subcell_cutPoints_global = nodes_n_cutpoints_global(active_nodes(:,:,icut),:,icut);
                
                [new_unfitted_cut_coord_iso,new_unfitted_cut_coord_global,new_x_unfitted_cut,new_subcell_cut_interior_connec_iso]...
                    = obj.computeInteriorSubcells(obj.fitted_mesh.connec(icell,:),subcell_cutPoints_iso,subcell_cutPoints_global);         
                
                new_nodes_containing_cell = repmat(icell,[size(new_unfitted_cut_coord_iso,1) 1]);
                
                k1 = k0 + size(new_unfitted_cut_coord_iso,1);
                obj.unfitted_cut_coord_iso(1+k0:k1,:) = new_unfitted_cut_coord_iso;
                obj.unfitted_cut_coord_global_raw(1+k0:k1,:) = new_unfitted_cut_coord_global;
                obj.nodes_containing_cell(1+k0:k1,:) = new_nodes_containing_cell;
                obj.x_unfitted_cut(1+k0:k1) = new_x_unfitted_cut;
                k0 = k1;
                
                new_subcell_containing_cell = repmat(icell,[size(new_subcell_cut_interior_connec_iso,1) 1]);
                m1 = m0+size(new_subcell_cut_interior_connec_iso,1);
                obj.unfitted_cut_connec_iso(1+m0:m1,:) = new_subcell_cut_interior_connec_iso;
                obj.subcell_containing_cell(1+m0:m1,:) = new_subcell_containing_cell;
                m0 = m1;
                
                c1 = c0 + size(new_subcell_cut_interior_connec_iso,1);
                obj.assignUnfittedCutCoordIsoPerCell(new_unfitted_cut_coord_iso,new_subcell_cut_interior_connec_iso,c0,c1);
                c0 = c1;
            end
            
            obj.cleanExtraAllocatedMemory(k1,m1,c1);
        end        
        function computeDelaunay_allocateMemory(obj)
            obj.unfitted_cut_coord_iso = zeros(length(obj.cut_cells)*obj.max_subcells*obj.nnodes_subcell,obj.fitted_mesh.ndim);
            obj.unfitted_cut_coord_global_raw = zeros(length(obj.cut_cells)*obj.max_subcells*obj.nnodes_subcell,obj.fitted_mesh.ndim);
            obj.x_unfitted_cut = zeros(length(obj.cut_cells)*obj.max_subcells*obj.nnodes_subcell,1);
            obj.unfitted_cut_connec_iso = zeros(length(obj.cut_cells)*obj.max_subcells,obj.nnodes_subcell);
            obj.unfitted_cut_coord_iso_per_cell = zeros(length(obj.cut_cells)*obj.max_subcells,obj.nnodes_subcell,obj.fitted_mesh.ndim);
            obj.unfitted_cut_connec_global = zeros(length(obj.cut_cells)*obj.max_subcells,obj.nnodes_subcell);
            obj.subcell_containing_cell = zeros(length(obj.cut_cells)*obj.max_subcells*obj.nnodes_subcell,1);
        end
        function subcell_connec = computeInteriorSubcellsConnectivities(obj,subcell_coord_iso,subcell_x_value)
            DT = delaunayTriangulation(subcell_coord_iso);
            subcell_connec = DT.ConnectivityList;
            is_interior = all(subcell_x_value(subcell_connec) <= 0,2);
            if obj.flag_debug
                if any(obj.cut_cells(obj.i_debug)==obj.debug_elems)
                    f = figure('visible','off');
                    
                    triplot(DT)
                    axis off
                    for m=1:size(subcell_coord_iso,1)
                        text(subcell_coord_iso(m,1),subcell_coord_iso(m,2),num2str(subcell_x_value(m)));
                    end
                    g=sprintf('Elem: %d ', obj.cut_cells(obj.i_debug));
                    title(g);
                    saveas(f,fullfile(strcat('Delaunay_',num2str(obj.cut_cells(obj.i_debug)),'.png')))
                end
                obj.i_debug=obj.i_debug+1;
            end
            subcell_connec = subcell_connec(is_interior,:);
        end
    end    
end

