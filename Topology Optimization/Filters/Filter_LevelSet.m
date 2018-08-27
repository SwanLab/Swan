classdef Filter_LevelSet < handle
    properties
        quadrature
        geometry
        quadrature_del
        interp_del
        shape_full
    end
    
    properties (Access = protected)
        max_subcells
        nnodes_subelem
        ndim
    end
    
    methods
        function preProcess(obj)
            obj.quadrature = Quadrature.set(obj.diffReacProb.geometry.type);
            obj.geometry= Geometry(obj.diffReacProb.mesh,'LINEAR');
            
            obj.getQuadratureDel;
            obj.quadrature_del.computeQuadrature('LINEAR');
            mesh_del = obj.getMeshDel;
            obj.getInterpolationDel(mesh_del);
            obj.interp_del.computeShapeDeriv(obj.quadrature_del.posgp)
            
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

            subcells_coord=zeros(length(cut_elem)*obj.max_subcells,obj.nnodes_subelem,obj.ndim);
            phi_cut=zeros(length(cut_elem)*obj.max_subcells,obj.nnodes_subelem);
            global_connec=zeros(length(cut_elem)*obj.max_subcells,1);
            
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
            if ~isempty(cut_elem)
                [subcells_coord,global_connec,phi_cut]=obj.computeDelaunay(x,cut_elem);
                dvolu_cut=obj.computeDvoluCut(subcells_coord);
                pos_gp_del_natural=obj.computePosGpDelaunayNatural(subcells_coord);
                obj.geometry.interpolation.computeShapeDeriv(pos_gp_del_natural');
                shape_all=obj.integrateCut(phi_cut, global_connec, dvolu_cut, shape_all);
            end
            M2=obj.rearrangeOutputRHS(shape_all);
        end
        
        % !!!!!!!!!!!!!!!!!! REMOVED M2=computeRHS_facet !!!!!!!!!!!!!!!!!!
        
        function S = computeFacetSurface(obj,x)
            M2 = obj.computeRHS_facet(x,ones(size(x)));
            S = sum(M2);
        end
    end
    
    methods (Abstract)
        getQuadratureDel(obj)
        getMeshDel(obj)
        getInterpolationDel(obj,mesh_del)
        computeRHS_facet(obj,x,F)
        findCutPoints_Iso(obj,x,cut_elem)
%         findCutPoints_Global(obj,x,cut_elem)
%         createFacet(obj)
        computeDvoluCut(elcrd)
%         mapping(elem_cutPoints_global,facets_connectivities,facet_deriv,dvolu)
    end
end