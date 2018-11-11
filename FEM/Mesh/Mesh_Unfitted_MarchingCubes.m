classdef Mesh_Unfitted_MarchingCubes < Mesh_Unfitted
    properties
    end
    methods
        function obj = computeCutMesh(obj)
            [nodes_n_cutpoints_iso,cut_to_elem_connec] = obj.findCutPoints_Iso;
            cases_connec=obj.findMarchingCubesCases;
            if obj.flag_debug
                for i_elem=1:length(obj.cut_cells)
                    if any(obj.cut_cells(i_elem)==obj.debug_elems)
                        elem_nodes=cut_to_elem_connec(i_elem,:)';
                        connec=cases_connec(:,:,i_elem);
                        f = figure('visible','off');
                        
                        x_coords=nodes_n_cutpoints_iso(elem_nodes(elem_nodes>0),1);
                        y_coords=nodes_n_cutpoints_iso(elem_nodes(elem_nodes>0),2);
                        x_values=[obj.x_fitted(obj.fitted_mesh.connec(obj.cut_cells(i_elem),:));zeros(length(elem_nodes)-size(obj.fitted_mesh.connec,2),1)];
                        triplot(connec(any(connec>0,2),:),x_coords,y_coords)
                        axis off
                        for m=1:size(x_coords,1)
                            text(x_coords(m),y_coords(m),num2str(x_values(m)));
                        end
                        g=sprintf('Elem: %d ', obj.cut_cells(i_elem));
                        title(g);
                        saveas(f,fullfile(strcat('Marching_',num2str(obj.cut_cells(i_elem)),'.png')))
                    end
                end
            end
            
            [elecoord_main,obj.x_unfitted_cut,obj.subcell_containing_cell]=obj.findSubElemCoord(cases_connec,cut_to_elem_connec,nodes_n_cutpoints_iso);
            
            for idime=1:obj.fitted_geom_interpolation.ndime
                obj.unfitted_cut_coord_iso_per_cell(:,:,idime)=elecoord_main{idime}(:,:);
            end
        end
        function   [cases_connec]=findMarchingCubesCases(obj)
            num=repmat([1:size(obj.fitted_mesh.connec,2)],[size(obj.cut_cells) 1]);
            cutcases=obj.x_fitted(obj.fitted_mesh.connec(obj.cut_cells,:))<0;
            negative_nodes_id=num.*cutcases;
            sum_negative_nodes_id=sum(negative_nodes_id,2);
            sum_negative_nodes=sum(cutcases,2);
            
            ind_cases=sub2ind(size(obj.fitted_geom_interpolation.selectcases),sum_negative_nodes_id,sum_negative_nodes);
            cases_list=obj.fitted_geom_interpolation.selectcases(ind_cases);
            
            cases_connec=obj.fitted_geom_interpolation.cases(:,:,cases_list);
        end
        function [sub_elem_coord,phi_cut,global_connec]=findSubElemCoord(obj,cases_connec,cut_to_elem_connec,cut_points)
            x_elem=[obj.x_fitted(obj.fitted_mesh.connec(obj.cut_cells,:)),zeros(length(obj.cut_cells),size(cut_to_elem_connec,2)-obj.fitted_geom_interpolation.nnode)];
            elem_list=1:length(obj.cut_cells);
            phi_cut=[];global_connec=[];
            sub_elem_coord=cell(obj.fitted_geom_interpolation.ndime,1);
            for i=1:size(cases_connec(:,:,1),1)
                coord_local=cell(obj.fitted_geom_interpolation.ndime,1);
                phi_local=[];
                for j=1:size(cases_connec(:,:,1),2)
                    connectivities=squeeze(cases_connec(i,j,:));
                    is_main_case = connectivities~=0;
                    if any(connectivities)
                        ind=sub2ind(size(cut_to_elem_connec),elem_list(is_main_case),connectivities(is_main_case)');
                        phi_local=[phi_local,x_elem(ind)'];
                        for idime=1:obj.fitted_geom_interpolation.ndime
                            coord_local{idime}=[coord_local{idime}(:,:),cut_points(cut_to_elem_connec(ind),idime)];
                        end
                    end
                end
                isnegative=~any(phi_local>0,2);
                for idime=1:obj.fitted_geom_interpolation.ndime
                    sub_elem_coord{idime}=[sub_elem_coord{idime}(:,:);coord_local{idime}(isnegative,:)];
                end                 
                phi_cut=[phi_cut;phi_local(isnegative,:)];
                local_connec=obj.cut_cells(is_main_case);
                global_connec=[global_connec;local_connec(isnegative)];
            end
        end
    end
end

