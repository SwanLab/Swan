classdef Filter_LevelSet < handle
    properties
        geometry % !! NEEDED??? !!
        quadrature
        quadrature_unfitted
        interpolation_unfitted
        unfitted_mesh
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
            obj.geometry= Geometry(obj.mesh,'LINEAR');
            
            obj.getQuadrature_Unfitted;
            obj.quadrature_unfitted.computeQuadrature('LINEAR');
            obj.createUnfittedMesh;
            obj.getInterpolation_Unfitted;
            obj.interpolation_unfitted.computeShapeDeriv(obj.quadrature_unfitted.posgp)
            
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
            shape = zeros(size(obj.mesh.connec,1),size(obj.mesh.connec,2));
            for igauss = 1:size(obj.geometry.interpolation.shape,2)
                shape = shape+obj.geometry.interpolation.shape(:,igauss)'.*obj.geometry.dvolu(:,igauss);
            end
        end
        
        function shape_cut = integrateCut(obj,containing_cell,dvolu_cut)
            dvolu_frac = sum(obj.geometry.dvolu,2)/obj.geometry.interpolation.dvolu;
            shape_cut = obj.geometry.interpolation.shape'.*dvolu_cut.*dvolu_frac(containing_cell);
        end
        
        function M2 = rearrangeOutputRHS(obj,shape_all)
            M2 = zeros(obj.npnod,1);
            for inode = 1:obj.nnode
                M2 = M2+accumarray(obj.mesh.connec(:,inode),shape_all(:,inode),[obj.npnod,1],@sum,0);
            end
        end
        
        %         function M2=computeRHS_OLD(obj,x)
        %             [full_elem,cut_elem]=obj.findCutElements(x,obj.connectivities);
        %             shape_all=obj.computeFullElements(full_elem);
        %             if ~isempty(cut_elem)
        %                 [cut_subcells_coord,global_elem_index_of_each_cut_elem,phi_values]=obj.computeDelaunay(x,cut_elem,obj.connectivities,obj.geometry.interpolation);
        %                 dvolu_cut=obj.computeDvoluCut(cut_subcells_coord);
        %                 posgp_iso=obj.computePosGpDelaunayIsoparametric(cut_subcells_coord);
        %                 obj.geometry.interpolation.computeShapeDeriv(posgp_iso');
        %                 shape_all=obj.integrateCut(phi_values, global_elem_index_of_each_cut_elem, dvolu_cut, shape_all);
        %             end
        %             M2=obj.rearrangeOutputRHS(shape_all);
        %         end
        
        function M2 = computeRHS(obj,x)
            obj.unfitted_mesh.computeCutMesh(x);
            obj.unfitted_mesh.computeDvoluCut;
            
            posgp_iso = obj.computePosGP(obj.unfitted_mesh.unfitted_cut_coord_iso_per_cell,obj.interpolation_unfitted);
            obj.geometry.interpolation.computeShapeDeriv(posgp_iso');
            
            shape_cut = obj.integrateCut(obj.unfitted_mesh.subcell_containing_cell,obj.unfitted_mesh.dvolu_cut);
            shape_all = obj.assembleShapeValues(shape_cut);
            
            M2=obj.rearrangeOutputRHS(shape_all);
        end
        
        function shape_all = assembleShapeValues(obj,shape_cut)
            shape_all = zeros(size(obj.mesh.connec,1),size(obj.mesh.connec,2));
            shape_all(obj.unfitted_mesh.full_cells,:) = obj.shape_full(obj.unfitted_mesh.full_cells,:);
            
            for i_subcell=1:size(shape_cut,2)
                shape_all(:,i_subcell)=shape_all(:,i_subcell)+accumarray(obj.unfitted_mesh.subcell_containing_cell,shape_cut(:,i_subcell),[obj.nelem,1],@sum,0);
            end
        end
        
        % !!!!!!!!!!!!!!!!!! REMOVED M2=computeRHS_facet !!!!!!!!!!!!!!!!!!
        
        function S = IntegrateFacet(obj,x)
            M2 = obj.computeRHS_facet(x,ones(size(x)));
            S = sum(M2);
        end
        
        function S = IntegrateInteriorCells(obj,x)
            M2 = obj.computeRHS(x);
            S = sum(M2);
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
        
        function posgp = computePosGP(subcell_coord,interpolation)
            posgp = zeros(size(subcell_coord,1),size(subcell_coord,3));
            for idime = 1:size(subcell_coord,3)
                posgp(:,idime) = subcell_coord(:,:,idime)*interpolation.shape;
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