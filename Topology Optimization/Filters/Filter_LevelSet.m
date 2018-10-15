classdef Filter_LevelSet < handle
    properties (Access = protected)
        quadrature_fitted
        quadrature_unfitted
        interpolation_fitted
        interpolation_unfitted
        unfitted_mesh
    end
    
    properties (Access = protected) % !! TO REMOVE !!
        max_subcells
        nnodes_subelem
        ndim
    end
    
    methods (Abstract)
        createUnfittedMesh(obj)
        setInterpolation_Unfitted(obj)
    end
    
    methods (Access = public)
        function preProcess(obj)
            obj.setQuadrature_Fitted;
            obj.setInterpolation_Fitted;
            
            obj.setQuadrature_Unfitted;
            obj.createUnfittedMesh;
            obj.setInterpolation_Unfitted;
            
            MSGID = 'MATLAB:delaunayTriangulation:DupPtsWarnId';
            warning('off', MSGID)
        end
    end
    
    methods (Access = protected)
        function M2 = rearrangeOutputRHS(obj,shape_all)
            M2 = zeros(obj.npnod,1);
            for inode = 1:obj.nnode
                M2 = M2 + accumarray(obj.mesh.connec(:,inode),shape_all(:,inode),[obj.npnod,1],@sum,0);
            end
        end
    end
    
    methods (Access = private)        
        function setInterpolation_Fitted(obj)
            obj.interpolation_fitted = Interpolation.create(obj.mesh,obj.quadrature_fitted.order);
        end
        
        function setQuadrature_Fitted(obj)
            obj.quadrature_fitted = Quadrature.set(obj.diffReacProb.geometry.type);
            obj.quadrature_fitted.computeQuadrature('LINEAR');
        end
        
        function setQuadrature_Unfitted(obj)
            obj.quadrature_unfitted = obj.getQuadrature_Unfitted;
            obj.quadrature_unfitted.computeQuadrature('LINEAR');
        end
    end
    
    methods (Static, Access = protected)
        function [full_elem,cut_elem] = findCutElements(x,connectivities) % !! TO REMOVE !!
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
end