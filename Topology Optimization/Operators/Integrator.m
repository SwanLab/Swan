classdef Integrator < handle
    properties (GetAccess = public, SetAccess = private)
        mesh_unfitted
        mesh_background
    end
    
%     properties (GetAccess = protected, SetAccess = private)
%         interpolation_unfitted
%         interpolation_background
%         
%         quadrature_unfitted
%         quadrature_background
%     end       
    
    methods (Access = public)
        function A = integrateUnfittedMesh(obj,F,mesh_unfitted)
            if exist('mesh_unfitted','var')
                obj.saveMeshes(mesh_unfitted);
            end
            A = obj.computeIntegral(F);
        end
    end
    
    methods (Static, Access = public)
        function obj = create(mesh)
            obj = IntegratorFactory.create(mesh);
            obj.saveMeshes(mesh);
        end
    end
    
    methods (Access = protected)
        function shapeValues = integrateCutCells(obj,F1)
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
        
        function M2 = rearrangeOutputRHS(obj,shapeValues_AllCells)
            interpolation = Interpolation.create(obj.mesh_background,'LINEAR');
            
            M2 = zeros(interpolation.npnod,1);
            for inode = 1:interpolation.nnode
                M2 = M2 + accumarray(obj.mesh_background.connec(:,inode),shapeValues_AllCells(:,inode),[interpolation.npnod,1],@sum,0);
            end
        end
        
        function saveMeshes(obj,unfitted)
            background = unfitted.mesh_background;
            obj.mesh_unfitted = unfitted;
            obj.mesh_background = background;
        end
    end
    
    methods (Static, Access = protected)
        function quadrature = computeQuadrature(geometryType)
            quadrature = Quadrature.set(geometryType);
            quadrature.computeQuadrature('LINEAR');
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

