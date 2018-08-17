classdef Filter_LevelSet < handle
    properties
        quadrature
        geometry
        quadrature_del
        interp_del
        shape_full
    end
    
    methods
        function preProcess(obj)
            obj.quadrature = Quadrature.set(obj.diffReacProb.geometry.type);
            mesh=obj.diffReacProb.mesh;
            obj.geometry= Geometry(mesh,'LINEAR');
            mesh_del=mesh;
            switch mesh.pdim
                case '2D'
                    mesh_del.geometryType='TRIANGLE';
                    obj.quadrature_del=Quadrature_Triangle;
                    obj.quadrature_del.computeQuadrature('LINEAR');
                    obj.interp_del=Triangle_Linear(mesh_del);
                    obj.interp_del.computeShapeDeriv(obj.quadrature_del.posgp)
                case '3D'
                    mesh_del.geometryType='TETRAHEDRA';
                    obj.quadrature_del=Quadrature_Tetrahedra;
                    obj.quadrature_del.computeQuadrature('LINEAR');
                    obj.interp_del=Tetrahedra(mesh_del);
                    obj.interp_del.computeShapeDeriv(obj.quadrature_del.posgp)
            end
            obj.initGeometry
            obj.shape_full=obj.integrateFull;
           
            MSGID = 'MATLAB:delaunayTriangulation:DupPtsWarnId';
            warning('off', MSGID)
        end
        
        function initGeometry(obj)
            obj.quadrature.computeQuadrature('LINEAR');
            obj.geometry.interpolation.computeShapeDeriv(obj.quadrature.posgp)
            obj.geometry.computeGeometry(obj.quadrature,obj.geometry.interpolation);
        end
        
        function [full_elem,cut_elem]=findCutElements(obj,x)
            phi_nodes=x(obj.connectivities);
            phi_case=sum((sign(phi_nodes)<0),2);
            
            full_elem = phi_case==size(obj.connectivities,2);
            null_elem = phi_case==0;
            indexes = (1:size(obj.connectivities,1))';
            cut_elem = indexes(~(full_elem+null_elem));
        end
        
        function shape=integrateFull(obj)
            shape=zeros(size(obj.connectivities,1),size(obj.connectivities,2));
            for igauss=1:size(obj.geometry.interpolation.shape,2)
                shape=shape+obj.geometry.interpolation.shape(:,igauss)'.*obj.geometry.dvolu(:,igauss);
            end
        end
        
        function [P,active_nodes]=findCutPoints_Iso(obj,x,cut_elem)
            switch obj.diffReacProb.mesh.pdim
                case '2D'
                    gamma_1=permute(x(obj.connectivities(cut_elem,:)),[2 3 1]);
                    gamma_2=permute([x(obj.connectivities(cut_elem,2:end)),x(obj.connectivities(cut_elem,1))],[2 3 1]);
                    P1=repmat(obj.geometry.interpolation.pos_nodes,[1 1 size(cut_elem)]);
                    P2=repmat([obj.geometry.interpolation.pos_nodes(2:end,:);obj.geometry.interpolation.pos_nodes(1,:)],[1 1 size(cut_elem)]);
                    P=P1+gamma_1.*(P2-P1)./(gamma_1-gamma_2);
                    active_nodes = sign(gamma_1.*gamma_2)<0;
                case '3D'
                    iteration_1=obj.geometry.interpolation.iteration(1,:);
                    iteration_2=obj.geometry.interpolation.iteration(2,:);
                    gamma_1=permute(x(obj.connectivities(cut_elem,iteration_1)),[2 3 1]);
                    gamma_2=permute(x(obj.connectivities(cut_elem,iteration_2)),[2 3 1]);
                    P1=repmat(obj.geometry.interpolation.pos_nodes(iteration_1,:),[1 1 size(cut_elem)]);
                    P2=repmat(obj.geometry.interpolation.pos_nodes(iteration_2,:),[1 1 size(cut_elem)]);
                    P=P1+gamma_1.*(P2-P1)./(gamma_1-gamma_2);
                    active_nodes = sign(gamma_1.*gamma_2)<0;
            end
        end
        
        function [P,active_nodes]=findCutPoints_Global(obj,x,cut_elem)
            switch obj.diffReacProb.mesh.pdim
                case '2D'
                    index1 = permute(obj.connectivities(cut_elem,:),[2 3 1]);
                    index2 = [permute(obj.connectivities(cut_elem,2:end),[2 3 1]);...
                        permute(obj.connectivities(cut_elem,1),[2 3 1])];
                    gamma_1=x(index1);
                    gamma_2=x(index2);
                    coord1 = obj.coordinates(:,1); coord2 = obj.coordinates(:,2);
                    P1=[coord1(index1) coord2(index1)];
                    P2=[coord1(index2) coord2(index2)];
                    P=P1+gamma_1.*(P2-P1)./(gamma_1-gamma_2);
                    active_nodes = sign(gamma_1.*gamma_2)<0;
                case '3D'
                    iteration_1=obj.geometry.interpolation.iteration(1,:);
                    iteration_2=obj.geometry.interpolation.iteration(2,:);
                    
                    index1 = permute(obj.connectivities(cut_elem,iteration_1),[2 3 1]);
                    index2 = permute(obj.connectivities(cut_elem,iteration_2),[2 3 1]);
                    gamma_1=x(index1);
                    gamma_2=x(index2);
                    coord1 = obj.coordinates(:,1); coord2 = obj.coordinates(:,2); coord3 = obj.coordinates(:,3);
                    P1=[coord1(index1) coord2(index1) coord3(index1)];
                    P2=[coord1(index2) coord2(index2) coord3(index2)];
                    P=P1+gamma_1.*(P2-P1)./(gamma_1-gamma_2);
                    active_nodes = sign(gamma_1.*gamma_2)<0;
            end
        end
        
        function A=computeDvoluCut(obj,elcrd)
            switch obj.diffReacProb.mesh.pdim
                case '2D'
                    x1 = elcrd(:,1,1); y1 = elcrd(:,1,2); x2 = elcrd(:,2,1); y2 = elcrd(:,2,2); x3 = elcrd(:,3,1); y3 = elcrd(:,3,2);
                    A = 0.5*abs((x2-x1).*(y3-y1)-(x3-x1).*(y2-y1));
                case '3D'
                    x1 = elcrd(:,1,1); y1 = elcrd(:,1,2); z1=elcrd(:,1,3);
                    x2 = elcrd(:,2,1); y2 = elcrd(:,2,2); z2=elcrd(:,2,3);
                    x3 = elcrd(:,3,1); y3 = elcrd(:,3,2); z3=elcrd(:,3,3);
                    x4 = elcrd(:,4,1); y4 = elcrd(:,4,2); z4=elcrd(:,4,3);
                    J=x1.*y3.*z2-x1.*y2.*z3+x2.*y1.*z3-x2.*y3.*z1-x3.*y1.*z2+x3.*y2.*z1+x1.*y2.*z4-x1.*y4.*z2-x2.*y1.*z4+x2.*y4.*z1+...
                        x4.*y1.*z2-x4.*y2.*z1-x1.*y3.*z4+x1.*y4.*z3+x3.*y1.*z4-x3.*y4.*z1-x4.*y1.*z3+x4.*y3.*z1+x2.*y3.*z4-x2.*y4.*z3...
                        -x3.*y2.*z4+x3.*y4.*z2+x4.*y2.*z3-x4.*y3.*z2;
                    A=J/6;
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
        
        function pos_gp_del_natural=computePosGpDelaunayNatural(obj,elcrd)
            pos_gp_del_natural=zeros(size(elcrd,1),size(elcrd,3));
            for idime=1:size(elcrd,3)
                pos_gp_del_natural(:,idime)=elcrd(:,:,idime)*obj.interp_del.shape;
            end
        end
        
        function shape_all=computeFullElements(obj,full_elem)
            shape_all=zeros(size(obj.connectivities,1),size(obj.connectivities,2));
            shape_all(full_elem,:)=obj.shape_full(full_elem,:);
        end
        
        function [subcells_coord,global_connec,phi_cut]=computeDelaunay(obj,x,cut_elem)
            [P,active_nodes]=obj.findCutPoints_Iso(x,cut_elem);
            switch obj.diffReacProb.mesh.pdim
                case '2D'
                    max_subcells = 6; nnodes_subelem = 3; ndim = 2;
                case '3D'
                    max_subcells = 20; nnodes_subelem = 4; ndim = 3;
            end
            subcells_coord=zeros(length(cut_elem)*max_subcells,nnodes_subelem,ndim);
            phi_cut=zeros(length(cut_elem)*max_subcells,nnodes_subelem);
            global_connec=zeros(length(cut_elem)*max_subcells,1);
            k = 0; m0 = 0;
            for ielem=1:length(cut_elem)
                del_coord = [obj.geometry.interpolation.pos_nodes;P(active_nodes(:,:,ielem),:,ielem)];
                del_x=[x(obj.connectivities(cut_elem(ielem),:));zeros(size(P(active_nodes(:,:,ielem)),1),1)]';
                DT=delaunayTriangulation(del_coord);
                del_connec=DT.ConnectivityList;
                new_subcells_coord = permute(del_coord,[3 1 2]);
                for idelaunay=1:size(del_connec,1)
                    k = k+1;
                    subcells_coord(k,:,:) = new_subcells_coord(:,del_connec(idelaunay,:),:);
                    phi_cut(k,:) = del_x(del_connec(idelaunay,:));
                end
                new_global_connec = repmat(cut_elem(ielem),[size(del_connec,1) 1]);
                m1 = m0+length(new_global_connec);
                global_connec(1+m0:m1,:)=repmat(cut_elem(ielem),[size(del_connec,1) 1]);
                m0 = m1;
            end
            if length(subcells_coord) > k
                subcells_coord(k+1:end,:,:) = [];
                phi_cut(k+1:end,:) = [];
                global_connec(m1+1:end) = [];
            end
        end
        
        function M2=rearrangeOutputRHS(obj,shape_all)
            M2=zeros(obj.npnod,1);
            for inode=1:obj.nnode
                p = obj.connectivities(:,inode);
                M2 = M2+accumarray(p,shape_all(:,inode),[obj.npnod,1],@sum,0);
            end
        end
        
        function M2=computeRHS(obj,x)
            [full_elem,cut_elem]=obj.findCutElements(x);
            shape_all=obj.computeFullElements(full_elem);
            % !!!!!!!!!!!!!!!!!!!!!!!!!! DELETE !!!!!!!!!!!!!!!!!!!!!!!!!!!
            %             shape_all(cut_elem,:)=0.5*obj.shape_full(cut_elem,:);
            
            %             [P,active_nodes]=obj.findCutPoints_Global(x,cut_elem);
            %             P = permute(P,[1 3 2]);
            %             active_nodes = permute(active_nodes,[1 3 2]);
            
            %             X = P(:,:,1); Y = P(:,:,2); Z = P(:,:,3);
            %             figure, plot3(X(active_nodes), Y(active_nodes), Z(active_nodes),'.')
            
            if ~isempty(cut_elem)
                [subcells_coord,global_connec,phi_cut]=obj.computeDelaunay(x,cut_elem);
                dvolu_cut=obj.computeDvoluCut(subcells_coord);
                pos_gp_del_natural=obj.computePosGpDelaunayNatural(subcells_coord);
                obj.geometry.interpolation.computeShapeDeriv(pos_gp_del_natural');
                shape_all=obj.integrateCut(phi_cut, global_connec, dvolu_cut, shape_all);
            end
            M2=obj.rearrangeOutputRHS(shape_all);
        end
        
        function M2=computeRHS_facet(obj,x,F)
            [interp_facet,quadrature_facet] = obj.createFacet;
            interp_element = Interpolation.create(obj.diffReacProb.mesh,obj.quadrature.order);
            
            shape_all = zeros(obj.nelem,obj.nnode);
            [~,cut_elem]=obj.findCutElements(x);
            
            [P_iso,active_nodes_iso]=obj.findCutPoints_Iso(x,cut_elem);
            [P_global,active_nodes_global]=obj.findCutPoints_Global(x,cut_elem);
            
            % !! VECTORITZAR: LOOPS PETITS, ELEMENTS DIRECTES !!
            %             figure, hold on
            for icut = 1:length(cut_elem)
                ielem = cut_elem(icut); inode_global = obj.connectivities(ielem,:);
                cutPoints_iso = P_iso(active_nodes_iso(:,:,icut),:,icut);
                cutPoints_global = P_global(active_nodes_global(:,:,icut),:,icut);
                
                connec_facets = obj.findFacetsConnectivities(cutPoints_iso,interp_element,x,inode_global);
                
                for i = 1:size(connec_facets,1)
                    for igaus = 1:quadrature_facet.ngaus
                        for idime = 1:interp_element.ndime
                            facet_posgp(igaus,idime) = interp_facet.shape(igaus,:)*cutPoints_iso(connec_facets(i,:),idime);
                        end
                    end
                    interp_element.computeShapeDeriv(facet_posgp');
                    facet_deriv(:,:) = interp_facet.deriv(:,:,:);
                    
                    % !! How mapping is done for 2D cases??? !!
                    t = [0; norm(diff(cutPoints_global(connec_facets(i,:),:)))];
                    dt_dxi = (facet_deriv'*t)/interp_facet.dvolu;
                    
                    f = (interp_element.shape*quadrature_facet.weigp')'*F(inode_global)/interp_facet.dvolu;
                    shape_all(ielem,:) = shape_all(ielem,:) + (interp_element.shape*(dt_dxi.*quadrature_facet.weigp')*f)';
                    
                    %                     plot(obj.coordinates(obj.connectivities(ielem,:),1),obj.coordinates(obj.connectivities(ielem,:),2),'.-b'); plot(obj.coordinates(obj.connectivities(ielem,[1 4]),1),obj.coordinates(obj.connectivities(ielem,[1 4]),2),'.-b');
                    %                     plot(cutPoints_global(connec_facets(i,:),1),cutPoints_global(connec_facets(i,:),2),'-xr');
                    %                     title('Cut Elements & Cut points in GLOBAL coordinates'), axis('equal')
                end
            end
            M2=obj.rearrangeOutputRHS(shape_all);
        end
        
        function [interp_facet,quadrature_facet] = createFacet(obj)
            switch obj.diffReacProb.mesh.pdim
                case '2D'
                    quadrature_facet = Quadrature.set('LINE');
                    interp_facet = Line_Linear;
                case '3D'
                    error('facets still NOT implemented for 3D meshes');
            end
            quadrature_facet.computeQuadrature(obj.quadrature.order);
            interp_facet.computeShapeDeriv(quadrature_facet.posgp);
        end
    end
    
    methods (Static)
        function connec_facets = findFacetsConnectivities(cutPoints_iso,interpolation,x,inode_global)
            if size(cutPoints_iso,1) == 2
                connec_facets = [1 2];
            elseif size(cutPoints_iso,1) == 4
                del_coord = [interpolation.pos_nodes;cutPoints_iso];
                DT=delaunayTriangulation(del_coord);
                del_connec=DT.ConnectivityList;
                
                node_positive_iso = find(x(inode_global)>0);
                
                for idel = 1:length(node_positive_iso)
                    [a, ~] = find(del_connec==node_positive_iso(idel));
                    connec_facets(idel,:) = del_connec(a(end),del_connec(a(end),:)~=node_positive_iso(idel))-interpolation.nnode;
                end
            else
                error('Case still not implemented.')
            end
        end
    end
end