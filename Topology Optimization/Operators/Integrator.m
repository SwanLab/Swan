classdef Integrator < handle
    properties (GetAccess = public, SetAccess = private)
        mesh_unfitted
        mesh_background
        
        %         interpolation_unfitted
        %         interpolation_background
        %
        %         quadrature_unfitted
        %         quadrature_background
    end
    
    methods (Access = public)
        function shapeValues = integrateMesh(mesh,F1)
            shapeValues = 'eis';
        end
        
        function M2 = integrateUnfittedMesh(obj,mesh_unfitted,mesh_background,F1)
            obj.saveMeshes(mesh_unfitted,mesh_background);
            type = obj.identifyUnfittedMeshType(mesh_unfitted);
            switch type
                case 'INTERIOR'
                    M2 = obj.integrateInteriorDomain(mesh_unfitted,mesh_background,F1);
                case 'BOUNDARY'
                    M2 = obj.integrateBoundaryDomain(mesh_unfitted,mesh_background,F1);
            end
        end
        
        function M2 = integrateInteriorDomain(obj,mesh_unfitted,mesh_background,F1)
            shapeValues_CutCells = obj.integrateCutCells(mesh_unfitted,mesh_background,F1);
            shapeValues_FullCells = obj.integrateFullCells(mesh_background,F1); % !! shapeValues_FullCells could saved in a property instead of being re-computed all the time !!
            shapeValues_All = obj.assembleShapeValues_Interior(shapeValues_CutCells,shapeValues_FullCells);
            M2 = obj.rearrangeOutputRHS(shapeValues_All);
        end
        
        function M2 = integrateBoundaryDomain(obj,mesh_unfitted,mesh_background,F1)
            shapeValues = obj.integrateCutCells(mesh_unfitted,mesh_background,F1);
            shapeValues = obj.assembleShapeValues_Boundary(shapeValues);
            M2 = obj.rearrangeOutputRHS(shapeValues);
        end
    end
    
    methods (Access = private)
        function shapeValues = integrateCutCells(obj,mesh_unfitted,mesh_background,F1)
            interpolation_background = Interpolation.create(mesh_background,'LINEAR');
            interpolation_unfitted = Interpolation.create(mesh_unfitted,'LINEAR');
            quadrature_unfitted = obj.computeQuadrature(mesh_unfitted.geometryType);
            
            posGP_iso_unfitted = obj.computePosGP(mesh_unfitted.coord_iso_per_cell,interpolation_unfitted,quadrature_unfitted);
            
            shapeValues = zeros(size(mesh_unfitted.connec,1),interpolation_background.nnode);
            for isubcell = 1:size(mesh_unfitted.connec,1) % !! VECTORIZE THIS LOOP !!
                icell = mesh_unfitted.cell_containing_subcell(isubcell);
                inode = mesh_background.connec(icell,:);
                
                interpolation_background.computeShapeDeriv(posGP_iso_unfitted(:,:,isubcell)');
                
                djacob = obj.mapping(mesh_unfitted.coord(mesh_unfitted.connec(isubcell,:),:),interpolation_unfitted.dvolu); % !! Could be done through Geometry class?? !!
                
                F0 = (interpolation_background.shape*quadrature_unfitted.weigp')'*F1(inode)/interpolation_unfitted.dvolu;
                shapeValues(isubcell,:) = shapeValues(isubcell,:) + (interpolation_background.shape*(djacob.*quadrature_unfitted.weigp')*F0)';
            end
        end
        
        function shapeValues_FullCells = integrateFullCells(obj,mesh,F1) % !! Only used by when integrating Unfitted_Interior !!
            % !! F1 should be evalutated in Integrator and integration at Interior Full
            % Cells should be allowed !!
            
            interpolation = Interpolation.create(mesh,'LINEAR');
            quadrature = obj.computeQuadrature(mesh.geometryType);
            interpolation.computeShapeDeriv(quadrature.posgp);
            geometry = Geometry(mesh,'LINEAR');
            geometry.computeGeometry(quadrature,interpolation);
            
            shapeValues_FullCells = zeros(size(mesh.connec));
            for igauss = 1:quadrature.ngaus
                shapeValues_FullCells = shapeValues_FullCells + interpolation.shape(:,igauss)'.*geometry.dvolu(:,igauss);
                %                 shapeValues_FullCells = shapeValues_FullCells + interpolation.shape(:,igauss)'.*F1.*geometry.dvolu(:,igauss);
            end
        end
        
        function M2 = rearrangeOutputRHS(obj,shapeValues_AllCells)
            interpolation = Interpolation.create(obj.mesh_background,'LINEAR');
            
            M2 = zeros(interpolation.npnod,1);
            for inode = 1:interpolation.nnode
                M2 = M2 + accumarray(obj.mesh_background.connec(:,inode),shapeValues_AllCells(:,inode),[interpolation.npnod,1],@sum,0);
            end
        end
        
        function shapeValues_AllCells = assembleShapeValues_Boundary(obj,shapeValues_CutCells) % !! Only used by when integrating Unfitted_Boundary !!
            interpolation = Interpolation.create(obj.mesh_background,'LINEAR');
            shapeValues_AllCells = zeros(size(obj.mesh_background.connec));
            
            for i_subcell = 1:size(shapeValues_CutCells,2)
                shapeValues_AllCells(:,i_subcell) = shapeValues_AllCells(:,i_subcell)+accumarray(obj.mesh_unfitted.cell_containing_subcell,shapeValues_CutCells(:,i_subcell),[interpolation.nelem,1],@sum,0);
            end
        end
        
        function shapeValues_AllCells = assembleShapeValues_Interior(obj,shapeValues_CutCells,shapeValues_FullCells)
            interpolation = Interpolation.create(obj.mesh_background,'LINEAR');
            shapeValues_AllCells = zeros(size(obj.mesh_background.connec));
            shapeValues_AllCells(obj.mesh_unfitted.background_full_cells,:) = shapeValues_FullCells(obj.mesh_unfitted.background_full_cells,:);
            
            for i_subcell = 1:size(shapeValues_CutCells,2)
                shapeValues_AllCells(:,i_subcell) = shapeValues_AllCells(:,i_subcell)+accumarray(obj.mesh_unfitted.cell_containing_subcell,shapeValues_CutCells(:,i_subcell),[interpolation.nelem,1],@sum,0);
            end
        end
        
        function saveMeshes(obj,unfitted,background)
            obj.mesh_unfitted = unfitted;
            obj.mesh_background = background;
        end
    end
    
    methods (Static, Access = private)
        function type = identifyUnfittedMeshType(mesh_unfitted)
            if contains(class(mesh_unfitted),'BOUNDARY','IgnoreCase',true)
                type = 'BOUNDARY';
            else
                type = 'INTERIOR';
            end
        end
        
        function quadrature = computeQuadrature(geometryType)
            quadrature = Quadrature.set(geometryType);
            quadrature.computeQuadrature('LINEAR');
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
        
        function djacob = mapping(points,dvolu)
            % !! PERFORM THROUGH GEOMETRY CLASS OR EXTRACT "mapping" CAPACITY FROM
            % GEOMETRY CLASS !!
            
            N_points = size(points,1);
            switch N_points
                case 2
                    v = diff(points);
                    L = norm(v);
                    djacob = L/dvolu;
                case 3
                    if size(points,2) == 2
                        points = [points, zeros(N_points,1)];
                    end
                    v1 = diff(points([1 2],:));
                    v2 = diff(points([1 3],:));
                    A = 0.5*norm(cross(v1,v2));
                    djacob = A/dvolu;
                case 4
                    v1 = diff(points([1 2],:));
                    v2 = diff(points([1 3],:));
                    v3 = diff(points([1 4],:));
                    V = (1/6)*det([v1;v2;v3]);
                    djacob = V/dvolu;
            end
        end
    end
end

