classdef Filter_LevelSet_Boundary < Filter_LevelSet
    methods (Access = public)
        function preProcess(obj)
            preProcess@Filter_LevelSet(obj);
        end
        
        function M2 = computeRHS(obj,F1)
            facet_posgp_iso = obj.computePosGP(obj.unfitted_mesh.coord_iso_per_cell,obj.interpolation_unfitted,obj.quadrature_unfitted);
            
            shapeValues_OLD = zeros(obj.nelem,obj.nnode);
            for ifacet = 1:size(obj.unfitted_mesh.connec,1) % !! VECTORIZE THIS LOOP !!
                icell = obj.unfitted_mesh.cell_containing_subcell(ifacet);
                inode = obj.mesh.connec(icell,:);
                
                obj.interpolation.computeShapeDeriv(facet_posgp_iso(:,:,ifacet)');
                
                djacob = obj.mapping(obj.unfitted_mesh.coord(obj.unfitted_mesh.connec(ifacet,:),:),obj.interpolation_unfitted.dvolu); % !! Could be done through Geometry class?? !!
                
                F0 = (obj.interpolation.shape*obj.quadrature_unfitted.weigp')'*F1(inode)/obj.interpolation_unfitted.dvolu;
                shapeValues_OLD(icell,:) = shapeValues_OLD(icell,:) + (obj.interpolation.shape*(djacob.*obj.quadrature_unfitted.weigp')*F0)';
            end
            
            shapeValues = obj.integrateFoverMesh(F1);
            shapeValues = obj.assembleShapeValues(shapeValues);
            M2 = obj.rearrangeOutputRHS(shapeValues);
        end
        
        %         function M2 = computeRHS(obj,F1)
        %             facet_posgp_iso = obj.computePosGP(obj.unfitted_mesh.coord_iso_per_cell,obj.interpolation_unfitted,obj.quadrature_unfitted);
        %
        %             shapeValues = zeros(obj.nelem,obj.nnode);
        %             for ifacet = 1:size(obj.unfitted_mesh.connec,1) % !! VECTORIZE THIS LOOP !!
        %                 icell = obj.unfitted_mesh.cell_containing_subcell(ifacet);
        %                 inode = obj.mesh.connec(icell,:);
        %
        %                 obj.interpolation.computeShapeDeriv(facet_posgp_iso(:,:,ifacet)');
        %
        %                 djacob = obj.mapping(obj.unfitted_mesh.coord(obj.unfitted_mesh.connec(ifacet,:),:),obj.interpolation_unfitted.dvolu); % !! Could be done through Geometry class?? !!
        %
        %                 F0 = (obj.interpolation.shape*obj.quadrature_unfitted.weigp')'*F1(inode)/obj.interpolation_unfitted.dvolu;
        %                 shapeValues(icell,:) = shapeValues(icell,:) + (obj.interpolation.shape*(djacob.*obj.quadrature_unfitted.weigp')*F0)';
        %             end
        %             M2 = obj.rearrangeOutputRHS(shapeValues);
        %         end
        
        function S = computeSurface(obj,x)
            obj.unfitted_mesh.computeMesh(x);
            M2 = obj.computeRHS(ones(size(x)));
            S = sum(M2);
        end
    end
    
    methods (Access = private)
        function shapeValues_AllCells = assembleShapeValues(obj,shapeValues_CutCells)
            shapeValues_AllCells = zeros(size(obj.mesh.connec,1),size(obj.mesh.connec,2));

            for i_subcell = 1:size(shapeValues_CutCells,2)
                shapeValues_AllCells(:,i_subcell) = shapeValues_AllCells(:,i_subcell)+accumarray(obj.unfitted_mesh.cell_containing_subcell,shapeValues_CutCells(:,i_subcell),[obj.nelem,1],@sum,0);
            end
        end
    end
end