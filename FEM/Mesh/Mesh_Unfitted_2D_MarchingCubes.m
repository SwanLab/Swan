classdef Mesh_Unfitted_2D_MarchingCubes < Mesh_Unfitted_2D & Mesh_Unfitted_MarchingCubes
    
    methods
        function obj = Mesh_Unfitted_2D_MarchingCubes(fitted_mesh,x_fitted,fitted_geom_interpolation)
            obj@Mesh_Unfitted_2D(fitted_mesh,x_fitted,fitted_geom_interpolation);
        end
        function [P,global_connec] = findCutPoints_Iso(obj)
            iteration_1=1:size(obj.fitted_mesh.connec,2);
            iteration_2=[2:size(obj.fitted_mesh.connec,2),1];
            gamma_1=obj.x_fitted(obj.fitted_mesh.connec(obj.cut_cells,iteration_1));
            gamma_2=obj.x_fitted(obj.fitted_mesh.connec(obj.cut_cells,iteration_2));
            active_nodes = sign(gamma_1.*gamma_2)<0;
            P=[];
            last=size(obj.fitted_mesh.connec,2);
            list_elem=1:length(obj.cut_cells);
            global_connec=repmat(1:size(obj.fitted_mesh.connec,2),[size(obj.cut_cells,1) 1]);
            for iedge=1:size(iteration_1,2)
                gamma_1active=gamma_1(active_nodes(:,iedge),iedge);
                gamma_2active=gamma_2(active_nodes(:,iedge),iedge);
                P1=repmat(obj.fitted_geom_interpolation.pos_nodes(iteration_1(iedge),:),[size(gamma_1active) 1]);
                P2=repmat(obj.fitted_geom_interpolation.pos_nodes(iteration_2(iedge),:),[size(gamma_2active) 1]);
                P=[P;P1+gamma_1active.*(P2-P1)./(gamma_1active-gamma_2active)];
                list_points=(last+1):(last+size(gamma_1active,1));
                last=list_points(end);
                connec=list_elem(active_nodes(:,iedge))';
                global_connec(connec,iedge+size(obj.fitted_mesh.connec,2))=list_points';
            end
            for j=1:length(iteration_1)-1
                for i=(size(obj.fitted_mesh.connec,2)+1):(size(global_connec,2)-1)
                    change=global_connec(:,i)==0;
                    global_connec(change,i)=global_connec(change,i+1);
                    global_connec(change,i+1)=zeros(sum(change),1);
                end
            end
            P=[obj.fitted_geom_interpolation.pos_nodes;P];           
        end
    end
end

