classdef Filter_LevelSet < handle
    properties
        geometry
        quadrature_fitted
        quadrature_unfitted
        interpolation_unfitted
        unfitted_mesh
    end
    
    properties (Access = protected) % !! TO REMOVE !!
        max_subcells
        nnodes_subelem
        ndim
    end
    
    methods
        function preProcess(obj)
            obj.quadrature_fitted = Quadrature.set(obj.diffReacProb.geometry.type);
            obj.quadrature_fitted.computeQuadrature('LINEAR');
            
            obj.getQuadrature_Unfitted;
            obj.quadrature_unfitted.computeQuadrature('LINEAR');
            
            obj.createUnfittedMesh;
            
            obj.getInterpolation_Unfitted;
            
            obj.computeGeometry;            
            
            MSGID = 'MATLAB:delaunayTriangulation:DupPtsWarnId';
            warning('off', MSGID)
        end
        
        function computeGeometry(obj)
            obj.geometry = Geometry(obj.mesh,'LINEAR');
            obj.geometry.interpolation.computeShapeDeriv(obj.quadrature_fitted.posgp);
            obj.geometry.computeGeometry(obj.quadrature_fitted,obj.geometry.interpolation);
        end        
        
        function M2 = rearrangeOutputRHS(obj,shape_all)
            M2 = zeros(obj.npnod,1);
            for inode = 1:obj.nnode
                M2 = M2 + accumarray(obj.mesh.connec(:,inode),shape_all(:,inode),[obj.npnod,1],@sum,0);
            end
        end
    end
    
    methods (Static)
        function [full_elem,cut_elem] = findCutElements(x,connectivities)
            phi_nodes = x(connectivities);
            phi_case = sum((sign(phi_nodes)<0),2);
            
            full_elem = phi_case==size(connectivities,2);
            null_elem = phi_case==0;
            indexes = (1:size(connectivities,1))';
            cut_elem = indexes(~(full_elem+null_elem));
        end
        
        function posgp = computePosGP(subcell_coord,interpolation,quadrature)
            interpolation.computeShapeDeriv(quadrature.posgp);
            posgp = zeros(quadrature.ngaus,size(subcell_coord,3),size(subcell_coord,1));
            for igaus = 1:quadrature.ngaus
                for idime = 1:size(subcell_coord,3)
                    posgp(igaus,idime,:) = subcell_coord(:,:,idime)*interpolation.shape(:,igaus);
                end
            end
        end
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
end