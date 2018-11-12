classdef Filter_LevelSet < handle
    properties
        quadrature
        geometry
        quadrature_del
        interp_del
        shape_full
        unfitted_mesh
        unfitted_mesh_algorithm
    end
    
    properties (Access = protected)
        max_subcells
        nnodes_subelem
        ndim
    end
    
    methods
        function obj = Filter_LevelSet(unfitted_mesh_algorithm)
            obj.unfitted_mesh_algorithm=unfitted_mesh_algorithm;
        end
        function preProcess(obj)
            obj.quadrature = Quadrature.set(obj.diffReacProb.geometry.type);
            obj.geometry= Geometry(obj.diffReacProb.mesh,'LINEAR');
            
            obj.getQuadratureDel;
            obj.quadrature_del.computeQuadrature('QUADRATIC');
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
        
        function shape = integrateFull(obj)
            shape = zeros(size(obj.connectivities,1),size(obj.connectivities,2));
            for igauss = 1:size(obj.geometry.interpolation.shape,2)
                shape = shape+obj.geometry.interpolation.shape(:,igauss)'.*obj.geometry.dvolu(:,igauss);
            end
        end
        
        function shape_cut = integrateCut(obj,containing_cell,dvolu_cut,pos_gp_del_natural)
            dvolu_frac = sum(obj.geometry.dvolu,2)/obj.geometry.interpolation.dvolu;            
            shape_cut=zeros(size(containing_cell,1),obj.geometry.interpolation.nnode);
            for igauss=1:obj.quadrature_del.ngaus
                obj.geometry.interpolation.computeShapeDeriv(pos_gp_del_natural(:,:,igauss)');
                shape_cut = shape_cut+obj.geometry.interpolation.shape'.*dvolu_cut.*dvolu_frac(containing_cell)*obj.quadrature_del.weigp(igauss);
            end
        end
        
        function pos_gp_cut_natural = computePosGpCutElemNatural(obj,subcell_coord)
            pos_gp_cut_natural = zeros(size(subcell_coord,1),size(subcell_coord,3),size(obj.interp_del.shape,2));
            for igauss=1:size(obj.interp_del.shape,2)
                for idime = 1:size(subcell_coord,3)
                    pos_gp_cut_natural(:,idime,igauss) = subcell_coord(:,:,idime)*obj.interp_del.shape(:,igauss);
                end
            end
        end
        
        function M2=rearrangeOutputRHS(obj,shape_all)
            M2=zeros(obj.npnod,1);
            for inode=1:obj.nnode
                p = obj.connectivities(:,inode);
                M2 = M2+accumarray(p,shape_all(:,inode),[obj.npnod,1],@sum,0);
            end
        end
               
        function M2 = computeRHS(obj,x)
            obj.setupUnfittedMesh(x);            
            shape_global=obj.initializeGlobalNodalShapeFunctionValues;
            if obj.checkIfMeshIsCut
                obj.unfitted_mesh.computeCutMesh;
                obj.unfitted_mesh.computeDvoluCut;
                pos_gp_cut_natural = obj.computePosGpCutElemNatural(obj.unfitted_mesh.unfitted_cut_coord_iso_per_cell);
                shape_cut = obj.integrateCut(obj.unfitted_mesh.subcell_containing_cell,obj.unfitted_mesh.dvolu_cut,pos_gp_cut_natural);
                shape_global = obj.assembleCutShapeValuesToGlobal(shape_cut,shape_global);
            end                
            M2=obj.rearrangeOutputRHS(shape_global);
        end
        function shape_all = initializeGlobalNodalShapeFunctionValues(obj)
            shape_all = zeros(size(obj.connectivities,1),size(obj.connectivities,2));
            shape_all(obj.unfitted_mesh.full_cells,:) = obj.shape_full(obj.unfitted_mesh.full_cells,:);
        end        
        function shape_all = assembleCutShapeValuesToGlobal(obj,shape_cut,shape_all)            
            for idelaunay=1:size(shape_cut,2)
                shape_all(:,idelaunay)=shape_all(:,idelaunay)+accumarray(obj.unfitted_mesh.subcell_containing_cell,shape_cut(:,idelaunay),[obj.nelem,1],@sum,0);
            end
        end

        % !!!!!!!!!!!!!!!!!! REMOVED M2=computeRHS_facet !!!!!!!!!!!!!!!!!!
        function is_cut=checkIfMeshIsCut(obj)
            is_cut=length(obj.unfitted_mesh.cut_cells) > 1;
        end
        function S = computeFacetSurface(obj,x)
            M2 = obj.computeRHS_facet(x,ones(size(x)));
            S = sum(M2);
        end
        
        function V = computeInteriorVolume(obj,x)
            M2 = obj.computeRHS(x);
            V = sum(M2);
        end
    end
    
    methods (Static)
        function [full_elem,cut_elem]=findCutElements(x,connectivities)
            phi_nodes=x(connectivities);
            phi_case=sum((sign(phi_nodes)<0),2);
            
            full_elem = phi_case==size(connectivities,2);
            null_elem = phi_case==0;
            indexes = (1:size(connectivities,1))';
            cut_elem = indexes(~(full_elem+null_elem));
        end
        
        
        %     methods (Abstract)
        %         getQuadratureDel(obj)
        %         getMeshDel(obj)
        %         getInterpolationDel(obj,mesh_del)
        %         computeRHS_facet(obj,x,F)
        %         findCutPoints_Iso(obj,x,cut_elem,interpolation)
        %         %         findCutPoints_Global(obj,x,cut_elem)
        %         %         createFacet(obj)
        %         computeDvoluCut(elcrd)
        %         %         mapping(elem_cutPoints_global,facets_connectivities,facet_deriv,dvolu)
        %     end
        function [F,aire] = faireF2(p,t,psi)
            np = size(p,2); nt = size(t,2);
            F = zeros(np,1);
            p1 = t(1,:); p2 = t(2,:); p3 = t(3,:);
            x1 = p(1,p1); y1 = p(2,p1); x2 = p(1,p2); y2 = p(2,p2); x3 = p(1,p3); y3 = p(2,p3);
            A = 0.5*abs((x2-x1).*(y3-y1)-(x3-x1).*(y2-y1));
            
            beta = (psi<0);
            beta = pdeintrp(p,t,beta);
            k = find(beta>0.5);
            F = F+accumarray(p1(k)',A(k)/3',[np,1],@sum,0);
            F = F+accumarray(p2(k)',A(k)/3',[np,1],@sum,0);
            F = F+accumarray(p3(k)',A(k)/3',[np,1],@sum,0);
            aire = sum(A(k));
            
            k = find(abs(beta-1/3)<0.01);
            p1 = t(1,k); p2 = t(2,k); p3 = t(3,k);
            psi1 = psi(p1)'; psi2 = psi(p2)'; psi3 = psi(p3)';
            [psis,is] = sort([psi1;psi2;psi3],1);
            is = is+3*ones(3,1)*[0:length(k)-1];
            pl = [p1;p2;p3]; ps = pl(is);
            x1 = p(1,ps(1,:)); y1 = p(2,ps(1,:)); x2 = p(1,ps(2,:)); y2 = p(2,ps(2,:)); x3 = p(1,ps(3,:)); y3 = p(2,ps(3,:));
            x12 = (psis(1,:).*x2-psis(2,:).*x1)./(psis(1,:)-psis(2,:));
            y12 = (psis(1,:).*y2-psis(2,:).*y1)./(psis(1,:)-psis(2,:));
            x13 = (psis(1,:).*x3-psis(3,:).*x1)./(psis(1,:)-psis(3,:));
            y13 = (psis(1,:).*y3-psis(3,:).*y1)./(psis(1,:)-psis(3,:));
            A = 0.5*abs(((x12-x1).*(y13-y1)-(x13-x1).*(y12-y1)));
            F = F+accumarray(ps(1,:)',((1+psis(2,:)./(psis(2,:)-psis(1,:))+psis(3,:)./(psis(3,:)-psis(1,:))).*A/3)',[np,1],@sum,0);
            F = F+accumarray(ps(2,:)',((psis(1,:)./(psis(1,:)-psis(2,:))).*A/3)',[np,1],@sum,0);
            F = F+accumarray(ps(3,:)',((psis(1,:)./(psis(1,:)-psis(3,:))).*A/3)',[np,1],@sum,0);
            aire = aire+sum(A);
            
            k = find(abs(beta-2/3)<0.01);
            p1 = t(1,k); p2 = t(2,k); p3 = t(3,k);
            psi1 = psi(p1)'; psi2 = psi(p2)'; psi3 = psi(p3)';
            [psis,is] = sort([psi1;psi2;psi3],1,'descend');
            is = is+3*ones(3,1)*[0:length(k)-1];
            pl = [p1;p2;p3]; ps = pl(is);
            x1 = p(1,ps(1,:)); y1 = p(2,ps(1,:)); x2 = p(1,ps(2,:)); y2 = p(2,ps(2,:)); x3 = p(1,ps(3,:)); y3 = p(2,ps(3,:));
            x12 = (psis(1,:).*x2-psis(2,:).*x1)./(psis(1,:)-psis(2,:));
            y12 = (psis(1,:).*y2-psis(2,:).*y1)./(psis(1,:)-psis(2,:));
            x13 = (psis(1,:).*x3-psis(3,:).*x1)./(psis(1,:)-psis(3,:));
            y13 = (psis(1,:).*y3-psis(3,:).*y1)./(psis(1,:)-psis(3,:));
            A = 0.5*abs(((x12-x1).*(y13-y1)-(x13-x1).*(y12-y1)));
            F = F-accumarray(ps(1,:)',((1+psis(2,:)./(psis(2,:)-psis(1,:))+psis(3,:)./(psis(3,:)-psis(1,:))).*A/3)',[np,1],@sum,0);
            F = F-accumarray(ps(2,:)',((psis(1,:)./(psis(1,:)-psis(2,:))).*A/3)',[np,1],@sum,0);
            F = F-accumarray(ps(3,:)',((psis(1,:)./(psis(1,:)-psis(3,:))).*A/3)',[np,1],@sum,0);
            aire = aire-sum(A);
        end
    end
end