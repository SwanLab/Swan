classdef Filter_LevelSet_Boundary < Filter_LevelSet
    methods
        %         function preProcess(obj)
        %             preProcess@Filter_LevelSet(obj);
        %         end
        
        function M2 = computeRHS(obj,x,F)
            obj.unfitted_mesh.computeMesh(x);
            obj.unfitted_mesh.computeGlobalConnectivities;
            
            [interpolation_facet,quadrature_facet] = obj.createFacet;
            interp_element = Interpolation.create(obj.mesh,obj.quadrature_fitted.order);
            
            facet_posgp_iso = obj.computePosGP(obj.unfitted_mesh.coord_iso_per_cell,interpolation_facet,quadrature_facet);
            
            shape_all = zeros(obj.nelem,obj.nnode);
            
            for ifacet = 1:size(obj.unfitted_mesh.connec,1)
                icell = obj.unfitted_mesh.cell_containing_subcell(ifacet);
                inode = obj.mesh.connec(icell,:);
                facet_posgp = facet_posgp_iso(:,:,ifacet);
                interp_element.computeShapeDeriv(facet_posgp');
                
                djacob = obj.mapping(obj.unfitted_mesh.coord(obj.unfitted_mesh.connec(ifacet,:),:),interpolation_facet.dvolu);
                
                f = (interp_element.shape*quadrature_facet.weigp')'*F(inode)/interpolation_facet.dvolu;
                shape_all(icell,:) = shape_all(icell,:) + (interp_element.shape*(djacob.*quadrature_facet.weigp')*f)';
            end
            M2 = obj.rearrangeOutputRHS(shape_all);
        end
        
        function S = computeSurface(obj,x)
            M2 = obj.computeRHS(x,ones(size(x)));
            S = sum(M2);
        end
    end
end