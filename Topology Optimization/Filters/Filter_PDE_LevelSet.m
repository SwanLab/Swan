classdef Filter_PDE_LevelSet < Filter_PDE
    properties
        quadrature
        geometry        
        quadrature_del
        interp_del
    end
    
    methods
        function obj = Filter_PDE_LevelSet(problemID,scale)
            obj@Filter_PDE(problemID,scale);
        end
       function preProcess(obj)
            preProcess@Filter_PDE(obj)
            obj.quadrature = Quadrature.set(obj.diffReacProb.element.geometry.type);
            obj.quadrature = Quadrature.set(obj.diffReacProb.geometry.type);
            mesh=obj.diffReacProb.mesh;
            obj.geometry= Geometry(mesh,'LINEAR');         
            obj.quadrature_del=Quadrature_Triangle;
            mesh_del=mesh;            
            mesh_del.geometryType='TRIANGLE';
            obj.quadrature_del.computeQuadrature('LINEAR');
            obj.interp_del=Triangle_Linear(mesh_del);
            obj.interp_del.computeShapeDeriv(obj.quadrature_del.posgp)
        end
        function rhs = integrate_L2_function_with_shape_function(obj,x)
            %rhs = obj.faireF2(obj.coordinates',obj.connectivities',x);
            rhs=obj.computeDelaunay(x);
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
        function [P,active_nodes]=findCrossPoints(obj,x,cut_elem)
            gamma_1=permute(x(obj.connectivities(cut_elem,:)),[2 3 1]);
            gamma_2=permute([x(obj.connectivities(cut_elem,2:end)),x(obj.connectivities(cut_elem,1))],[2 3 1]);
            P1=repmat(obj.geometry.interpolation.pos_nodes,[1 1 size(cut_elem)]);
            P2=repmat([obj.geometry.interpolation.pos_nodes(2:end,:);obj.geometry.interpolation.pos_nodes(1,:)],[1 1 size(cut_elem)]);
            P=P1+gamma_1.*(P2-P1)./(gamma_1-gamma_2);
            active_nodes = sign(gamma_1) ~= sign(gamma_2);
        end
        function A=computeDvoluCut(obj,x_coord,y_coord,z_coord)
            switch obj.diffReacProb.mesh.pdim
                case '2D'
                    x1 = x_coord(:,1); y1 = y_coord(:,1); x2 = x_coord(:,2); y2 = y_coord(:,2); x3 = x_coord(:,3); y3 = y_coord(:,3);
                    A = 0.5*abs((x2-x1).*(y3-y1)-(x3-x1).*(y2-y1));
                case '3D'
            end
        end
        function shape_all=integrateCut(obj,phi_cut, global_connec, A, shape_all)
            notcompute = phi_cut > 0;
            isnegative = ~any(notcompute');
            dvolu=sum(obj.geometry.dvolu,2)/obj.geometry.interpolation.dvolu;
            v=isnegative'.*(obj.geometry.interpolation.shape'.*A.*dvolu(global_connec));
            for idelaunay=1:size(v,2)
                shape_all(:,idelaunay)=shape_all(:,idelaunay)+accumarray(global_connec,v(:,idelaunay),[obj.nelem,1],@sum,0);
            end
        end
        
        function M2=computeDelaunay(obj,x)
            obj.initGeometry
            [full_elem,cut_elem]=obj.findCases(x);
            shape_all=obj.integrateFull(full_elem);
            if ~isempty(cut_elem)
                [P,active_nodes]=obj.findCrossPoints(x,cut_elem);
                x_coord=[];y_coord=[];z_coord=[];phi_cut=[];global_connec=[];
                for ielem=1:length(cut_elem)
                    delaunay_coord = [obj.geometry.interpolation.pos_nodes;P(active_nodes(:,:,ielem),:,ielem)];
                    delaunay_x=[x(obj.connectivities(cut_elem(ielem),:));zeros(size(P(active_nodes(:,:,ielem)),1),1)];
                    delaunay_connec=delaunay(delaunay_coord(:,1),delaunay_coord(:,2));                    
                    for idelaunay=1:size(delaunay_connec,1)                        
                        x_coord = [x_coord; delaunay_coord(delaunay_connec(idelaunay,:),1)'];
                        y_coord = [y_coord; delaunay_coord(delaunay_connec(idelaunay,:),2)'];
                        phi_cut =[phi_cut;delaunay_x(delaunay_connec(idelaunay,:))'];
                    end
                    global_connec=[global_connec;repmat(cut_elem(ielem),[size(delaunay_connec,1) 1])];                    
                end
                dvolu_cut=obj.computeDvoluCut(x_coord,y_coord,z_coord);
   
                pos_gp_del_natural=[x_coord*obj.interp_del.shape,y_coord*obj.interp_del.shape];
                obj.geometry.interpolation.computeShapeDeriv(pos_gp_del_natural');
                
                shape_all=obj.integrateCut(phi_cut, global_connec, dvolu_cut, shape_all);                
            end            
            M2=zeros(obj.npnod,1);
            for inode=1:obj.nnode
                p = obj.connectivities(:,inode);
                M2 = M2+accumarray(p,shape_all(:,inode),[obj.npnod,1],@sum,0);
            end
        end
    end
end