classdef Filter_PDE_LevelSet < Filter_PDE
    properties
        quadrature
        geometry        
        quadrature_del
    end
    
    methods
       function preProcess(obj,physicalProblem)
            preProcess@Filter_PDE(obj,physicalProblem)
            obj.quadrature = Quadrature.set(physicalProblem.element.geometry.type);
            mesh=physicalProblem.mesh;
            obj.geometry= Geometry(mesh,'LINEAR');
            obj.quadrature_del=Quadrature_Triangle;
        end
        function rhs = integrate_L2_function_with_shape_function(obj,x)
            %rhs = obj.faireF2(obj.coordinates',obj.connectivities',x);
            rhs=obj.computeDelaunay(x);
        end
%         function x_gp = getP0fromP1(obj,x)
%             if norm(x) == norm(obj.x)
%                 x_gp=obj.x_reg;
%             else
%                                 %M2=obj.faireF2(obj.coordinates',obj.connectivities',x);
%                               % x_gp = obj.P_operator*M2;
%                 obj.x=x;
%              
%                 M2=obj.computeDelaunay(x);
%                 x_gp = obj.P_operator*M2;
%                 obj.x_reg=x_gp;
%             end
%         end
%         function x_reg = getP1fromP0(obj,x)
%             gauss_sum=0;
%             for igauss=1:size(obj.M0,2)
%                 if size(x,2)==1
%                     gauss_sum=gauss_sum+obj.M0{igauss}*x;
%                 else
%                     gauss_sum=gauss_sum+obj.M0{igauss}*x(:,igauss);
%                 end
%             end
%             x_reg = obj.P_operator'*gauss_sum;
%         end
        function M2=computeDelaunay(obj,x)
            coord=obj.coordinates;
            connec=obj.connectivities;
            obj.quadrature.computeQuadrature('LINEAR');
            obj.geometry.interpolation.computeShapeDeriv(obj.quadrature.posgp)
            obj.geometry.computeGeometry(obj.quadrature,obj.geometry.interpolation);
         
            phi_nodes=x(connec);
            phi_case=sum((sign(phi_nodes)<0),2);
            
            shape_all=zeros(size(connec,1),size(connec,2));
            full_elem = phi_case==size(connec,2);
            null_elem = phi_case==0;
            indexes = (1:size(connec,1))';
            delaunay_elem = indexes(~(full_elem+null_elem));
            
            %shape_all(full_elem,:)=geometry.dvolu(full_elem,:);
            for igauss=1:size(obj.geometry.interpolation.shape,2)
                shape_all(full_elem,:)=shape_all(full_elem,:)+obj.geometry.interpolation.shape(:,igauss)'.*obj.geometry.dvolu(full_elem,igauss);
            end
            if ~isempty(delaunay_elem)
                for ielem=1:length(delaunay_elem)
                    %interpolacion y puntos de corte
                    gamma_1=x(connec(delaunay_elem(ielem),:));
                    gamma_2=x(connec(delaunay_elem(ielem),2:end));
                    P1=obj.geometry.interpolation.pos_nodes;
                    P2=obj.geometry.interpolation.pos_nodes(2:end,:);
                    gamma_2=[gamma_2;x(connec(delaunay_elem(ielem),1))];
                    P2=[P2;obj.geometry.interpolation.pos_nodes(1,:)];
                    P=P1+gamma_1.*(P2-P1)./(gamma_1-gamma_2);
                    unactive_node = sign(gamma_1) == sign(gamma_2);
                    P(unactive_node,:)=[];
                    %calculo nuevos coord y connec
                    delaunay_coord = [obj.geometry.interpolation.pos_nodes;P];
                    delaunay_x=[x(connec(delaunay_elem(ielem),:));zeros(size(P,1),1)];
                    delaunay_connec=delaunay(delaunay_coord(:,1),delaunay_coord(:,2));
                    notcompute = delaunay_x(:,:)>0;
                    isnegative=~any((notcompute(delaunay_connec))');
                    
                    obj.geometry.interpolation.computeShapeDeriv(delaunay_coord')
                    
                    mesh_del.coord=obj.geometry.interpolation.shape'*coord(connec(delaunay_elem(ielem),:)',:);
                    
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
                        shape_all(delaunay_elem(ielem),:)=shape_all(delaunay_elem(ielem),:)+v;
                    end
                    
%                                         figure(3)
%                                         hold on
%                                         x_global=coord(connec(delaunay_elem(ielem),:)',1);
%                                         y_global=coord(connec(delaunay_elem(ielem),:)',2);
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
%                                         text(x_global,y_global,num2str(x(connec(delaunay_elem(ielem),:)')))
%                                         %title(num2str(phi(delaunay_elem(ielem))))
%                                         close (3)
                end
            end
            
            M2=zeros(obj.npnod,1);
            for inode=1:obj.nnode
                p = connec(:,inode);
                M2 = M2+accumarray(p,shape_all(:,inode),[obj.npnod,1],@sum,0);
            end
        end
    end
end