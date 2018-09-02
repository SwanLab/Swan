classdef Mesh_Unfitted < Mesh
    properties
        fitted_mesh
        fitted_geom_interpolation
        
        x_fitted
        x_unfitted
        x_unfitted_cut
        
        full_cells
        empty_cells
        cut_cells
        
        coord_cut
        container_cell
        dvolu_cut
        
        max_subcells
        nnodes_subcell
    end
    
    methods
        function obj = Mesh_Unfitted(fitted_mesh,x_fitted,fitted_geom_interpolation)
            obj.fitted_mesh = fitted_mesh;
            obj.fitted_geom_interpolation = fitted_geom_interpolation;
            obj.x_fitted = x_fitted;
        end
        
        function computeCutMesh(obj)
            obj.findCutCells;
            obj.computeCutDelaunay;
        end
        
        function obj = computeCutDelaunay(obj)
            [nodes_n_cutpoints,active_nodes] = obj.findCutPoints_Iso;
            
            obj.computeDelaunay_allocateMemory;
            
            k = 0; m0 = 0; m1 = 0;
            for icut = 1:length(obj.cut_cells)
                icell = obj.cut_cells(icut);
                
                del_coord = [obj.fitted_geom_interpolation.pos_nodes; nodes_n_cutpoints(active_nodes(:,:,icut),:,icut)];
                del_x_value = [obj.x_fitted(obj.fitted_mesh.connec(icell,:)); zeros(size(nodes_n_cutpoints(active_nodes(:,:,icut)),1),1)]';
                
                DT = delaunayTriangulation(del_coord);
                del_connec = DT.ConnectivityList;
                
                new_coord_cut = permute(del_coord,[3 1 2]);
                for idelaunay = 1:size(del_connec,1)
                    k = k+1;
                    obj.coord_cut(k,:,:) = new_coord_cut(:,del_connec(idelaunay,:),:);
                    obj.x_unfitted_cut(k,:) = del_x_value(del_connec(idelaunay,:));
                end
                new_container_cell = repmat(icell,[size(del_connec,1) 1]);
                m1 = m0+length(new_container_cell);
                obj.container_cell(1+m0:m1,:) = repmat(icell,[size(del_connec,1) 1]);
                m0 = m1;
            end
            
            obj.cleanExtraAllocatedMemory(k,m1);
        end
        
        function computeDelaunay_allocateMemory(obj)
            obj.coord_cut = zeros(length(obj.cut_cells)*obj.max_subcells,obj.nnodes_subcell,obj.fitted_mesh.ndim);
            obj.x_unfitted_cut = zeros(length(obj.cut_cells)*obj.max_subcells,obj.nnodes_subcell);
            obj.container_cell = zeros(length(obj.cut_cells)*obj.max_subcells,1);
        end
        
        function cleanExtraAllocatedMemory(obj,k,m)
            if length(obj.coord_cut) > k
                obj.coord_cut(k+1:end,:,:) = [];
                obj.x_unfitted_cut(k+1:end,:) = [];
                obj.container_cell(m+1:end) = [];
            end
        end
        
        function findCutCells(obj)
            phi_nodes = obj.x_fitted(obj.fitted_mesh.connec);
            phi_case = sum((sign(phi_nodes)<0),2);
            
            obj.full_cells = phi_case == size(obj.fitted_mesh.connec,2);
            obj.empty_cells = phi_case == 0;
            indexes = (1:size(obj.fitted_mesh.connec,1))';
            obj.cut_cells = indexes(~(obj.full_cells | obj.empty_cells));
        end
    end
end

