classdef Filter_P1_LevelSet < Filter_P1
    properties
        quadrature
        geometry        
        quadrature_del
    end
    methods    
        function obj = Filter_P1_LevelSet(problemID,scale)
            obj@Filter_P1(problemID,scale);
        end
        function preProcess(obj)
            preProcess@Filter_P1(obj)
            obj.quadrature = Quadrature.set(obj.diffReacProb.geometry.type);
            mesh=obj.diffReacProb.mesh;
            obj.geometry= Geometry(mesh,'LINEAR');         
            obj.quadrature_del=Quadrature_Triangle;
        end
        function x_gp = getP0fromP1(obj,x)
            if norm(x) == norm(obj.x)
                x_gp=obj.x_reg;
            else
                %M2=obj.faireF2(obj.coordinates',obj.connectivities',x);
                % x_gp = obj.P_operator*M2;
                obj.x=x;
                
                M2=obj.computeDelaunay(x);
                x_gp = obj.P_operator*M2;
                obj.x_reg=x_gp;
            end
        end
        function initGeometry(obj)
            obj.quadrature.computeQuadrature('LINEAR');
            obj.geometry.interpolation.computeShapeDeriv(obj.quadrature.posgp)
            obj.geometry.computeGeometry(obj.quadrature,obj.geometry.interpolation);
        end
        function [full_elem,cut_elem]=findCases(obj,x)
            phi_nodes=x(obj.connectivities);
            phi_case=sum((sign(phi_nodes)<0),2);     
            
            full_elem = phi_case==size(obj.connectivities,2);
            null_elem = phi_case==0;
            indexes = (1:size(obj.connectivities,1))';
            cut_elem = indexes(~(full_elem+null_elem));
        end
        function shape_all=integrateFull(obj,full_elem)
            shape_all=zeros(size(obj.connectivities,1),size(obj.connectivities,2));
            for igauss=1:size(obj.geometry.interpolation.shape,2)
                shape_all(full_elem,:)=shape_all(full_elem,:)+obj.geometry.interpolation.shape(:,igauss)'.*obj.geometry.dvolu(full_elem,igauss);
            end
        end
        function P=findCrossPoints(obj,x,cut_elem)
            gamma_1=x(obj.connectivities(cut_elem,:));
            gamma_2=x(obj.connectivities(cut_elem,2:end));
            P1=obj.geometry.interpolation.pos_nodes;
            P2=obj.geometry.interpolation.pos_nodes(2:end,:);
            gamma_2=[gamma_2;x(obj.connectivities(cut_elem,1))];
            P2=[P2;obj.geometry.interpolation.pos_nodes(1,:)];
            P=P1+gamma_1.*(P2-P1)./(gamma_1-gamma_2);
            unactive_node = sign(gamma_1) == sign(gamma_2);
            P(unactive_node,:)=[];
        end
        function M2=computeDelaunay(obj,x)
            obj.initGeometry
            [full_elem,cut_elem]=obj.findCases(x);
            shape_all=obj.integrateFull(full_elem);
           
            
            %shape_all(full_elem,:)=geometry.dvolu(full_elem,:);
            
            if ~isempty(cut_elem)
                for ielem=1:length(cut_elem)
                    %puntos de corte
                    P=obj.findCrossPoints(x,cut_elem(ielem));
                    %calculo nuevos coord y connec
                    delaunay_coord = [obj.geometry.interpolation.pos_nodes;P];
                    delaunay_x=[x(obj.connectivities(cut_elem(ielem),:));zeros(size(P,1),1)];
                    delaunay_connec=delaunay(delaunay_coord(:,1),delaunay_coord(:,2));
                    notcompute = delaunay_x(:,:)>0;
                    isnegative=~any((notcompute(delaunay_connec))');
                    
                    obj.geometry.interpolation.computeShapeDeriv(delaunay_coord')
                    
                    mesh_del.coord=obj.geometry.interpolation.shape'*obj.coordinates(obj.connectivities(cut_elem(ielem),:)',:);
                    
                    mesh_del.connec=delaunay_connec;
                    mesh_del.geometryType='TRIANGLE';
                    geometry_del=Geometry(mesh_del,'LINEAR');
                    obj.quadrature_del.computeQuadrature('LINEAR');
                    interp_del=Triangle_Linear(mesh_del);
                    interp_del.computeShapeDeriv(obj.quadrature_del.posgp)
                    geometry_del.computeGeometry(obj.quadrature_del,interp_del);               
                    
                    %integracion puntos de gauss negros
                    for idelaunay=1:size(delaunay_connec,1)
                        
                        pos_gp_del_natural=interp_del.shape'*delaunay_coord(delaunay_connec(idelaunay,:)',:);
                        obj.geometry.interpolation.computeShapeDeriv(pos_gp_del_natural');
                        v=isnegative(idelaunay)*obj.geometry.interpolation.shape(:,1)'.*geometry_del.dvolu(idelaunay);
                        shape_all(cut_elem(ielem),:)=shape_all(cut_elem(ielem),:)+v;
                    end
                    

                end
            end
            
            M2=zeros(obj.npnod,1);
            for inode=1:obj.nnode
                p = obj.connectivities(:,inode);
                M2 = M2+accumarray(p,shape_all(:,inode),[obj.npnod,1],@sum,0);
            end
        end
    end
end

%                                         figure(3)
%                                         hold on
%                                         x_global=coord(connec(cut_elem(ielem),:)',1);
%                                         y_global=coord(connec(cut_elem(ielem),:)',2);
%                                         for i=1:size(x_global,1)
%                                             vec=[i;i+1];
%                                             if i+1 > size(x_global,1)
%                                                 vec(2)=1;
%                                             end
%                                             plot(x_global(vec,:),y_global(vec,:),'b')
%                                         end
%                     
%                                         for idel=1:size(delaunay_connec,1)
%                                             x_delaunay=mesh_del.coord(delaunay_connec(idel,:)',1);
%                                             y_delaunay=mesh_del.coord(delaunay_connec(idel,:)',2);
%                                             for i=1:size(x_delaunay,1)
%                                                 vec=[i;i+1];
%                                                 if i+1 > size(x_delaunay,1)
%                                                     vec(2)=1;
%                                                 end
%                                                 hold on
%                                                 plot(x_delaunay(vec,:),y_delaunay(vec,:),'-r')
%                     
%                                                 drawnow
%                                             end
%                                         end
%                                         text(x_global,y_global,num2str(x(connec(cut_elem(ielem),:)')))
%                                         %title(num2str(phi(cut_elem(ielem))))
%                                         close (3)