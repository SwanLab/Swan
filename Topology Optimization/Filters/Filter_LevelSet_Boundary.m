classdef Filter_LevelSet_Boundary < Filter_LevelSet
    methods (Access = public)
        function preProcess(obj)
            preProcess@Filter_LevelSet(obj);
        end
        
        function M2 = computeRHS(obj,F)
            facet_posgp_iso = obj.computePosGP(obj.unfitted_mesh.coord_iso_per_cell,obj.interpolation_unfitted,obj.quadrature_unfitted);
            
            shapeValues = zeros(obj.nelem,obj.nnode);
            for ifacet = 1:size(obj.unfitted_mesh.connec,1) % !! VECTORIZE THIS LOOP !!
                icell = obj.unfitted_mesh.cell_containing_subcell(ifacet);
                inode = obj.mesh.connec(icell,:);
                
                obj.interpolation.computeShapeDeriv(facet_posgp_iso(:,:,ifacet)');
                
                djacob = obj.mapping(obj.unfitted_mesh.coord(obj.unfitted_mesh.connec(ifacet,:),:),obj.interpolation_unfitted.dvolu); % !! Could be done through Geometry class?? !!
                
                f = (obj.interpolation.shape*obj.quadrature_unfitted.weigp')'*F(inode)/obj.interpolation_unfitted.dvolu;
                shapeValues(icell,:) = shapeValues(icell,:) + (obj.interpolation.shape*(djacob.*obj.quadrature_unfitted.weigp')*f)';
            end
            M2 = obj.rearrangeOutputRHS(shapeValues);
        end
        
        function S = computeSurface(obj,x)
            obj.unfitted_mesh.computeMesh(x);
            obj.unfitted_mesh.computeGlobalConnectivities;
            M2 = obj.computeRHS(ones(size(x)));
            S = sum(M2);
        end
    end
end