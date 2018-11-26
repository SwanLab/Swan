classdef Integrator_Boundary < Integrator % !! Rename to: Integrator_Unfitted !!
    properties (GetAccess = public, SetAccess = private)
        mesh_unfitted
        mesh_background
    end
    
    properties (Access = private)
%         mapper
    end

    %     properties (Access = private)
    %         interpolation_unfitted
    %         interpolation_background
    %
    %         quadrature_unfitted
    %         quadrature_background
    %     end
    
    
    methods (Access = public)
        function obj = Integrator_Boundary(mesh)
            obj.saveMesh(mesh);
        end
    end

    methods (Access = private)
        function shapeValues = integrate(obj,F1)
            interpolation_background = Interpolation.create(obj.mesh_background,'LINEAR');
            interpolation_unfitted = Interpolation.create(obj.mesh_unfitted,'LINEAR');
            quadrature_unfitted = obj.computeQuadrature(obj.mesh_unfitted.geometryType);
            
            posGP_iso_unfitted = obj.computePosGP(obj.mesh_unfitted.coord_iso_per_cell,interpolation_unfitted,quadrature_unfitted);
            
            shapeValues = zeros(size(obj.mesh_unfitted.connec,1),interpolation_background.nnode);
            for isubcell = 1:size(obj.mesh_unfitted.connec,1) % !! VECTORIZE THIS LOOP !!
                icell = obj.mesh_unfitted.cell_containing_subcell(isubcell);
                inode = obj.mesh_background.connec(icell,:);
                
                interpolation_background.computeShapeDeriv(posGP_iso_unfitted(:,:,isubcell)');
                
                djacob = obj.mapping(obj.mesh_unfitted.coord(obj.mesh_unfitted.connec(isubcell,:),:),interpolation_unfitted.dvolu); % !! Could be done through Geometry class?? !!
                
                F0 = (interpolation_background.shape*quadrature_unfitted.weigp')'*F1(inode)/interpolation_unfitted.dvolu;
                shapeValues(isubcell,:) = shapeValues(isubcell,:) + (interpolation_background.shape*(djacob.*quadrature_unfitted.weigp')*F0)';
            end
        end
    end
    
    methods (Static, Access = private)
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

