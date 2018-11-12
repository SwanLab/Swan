classdef Integrator
    properties (GetAccess = public, SetAccess = private)
    end
    
    methods (Access = public)
        function shapeValues = integrateMesh(mesh,F1)
            shapeValues = 'eis';
        end
        
        
        function shapeValues = integrateUnfittedMesh(obj,mesh_unfitted,mesh_background,F1)
            interpolation_background = Interpolation.create(mesh_background,'LINEAR');
            interpolation_unfitted = Interpolation.create(mesh_unfitted,'LINEAR');
            quadrature_unfitted = Quadrature.set(mesh_unfitted.geometryType);
            quadrature_unfitted.computeQuadrature('LINEAR');
            
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
    end
    
    methods (Static, Access = private)
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

