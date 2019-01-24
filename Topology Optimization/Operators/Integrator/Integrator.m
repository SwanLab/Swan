classdef Integrator < handle
    
    properties (GetAccess = public, SetAccess = protected)
        meshUnfitted
        meshBackground
    end
    
    
    methods (Access = protected, Abstract)
        
        computeIntegral(obj)
        
    end
    
    methods (Access = public)
        
        function A = integrateUnfittedMesh(obj,F,meshUnfitted)
            if exist('meshUnfitted','var')
                obj.updateMeshes(meshUnfitted);
            end
            A = obj.computeIntegral(F);
        end
        
    end
    
    methods (Static, Access = public)
        
        function obj = create(mesh)
            obj = IntegratorFactory.create(mesh);
            obj.updateMeshes(mesh);
        end
        
    end
    
    methods (Access = protected)
        
        function shapeValues = integrateCutCells(obj,F1)
            interpolation_background = Interpolation.create(obj.meshBackground,'LINEAR');
            interpolation_unfitted = Interpolation.create(obj.meshUnfitted,'LINEAR');
            quadrature_unfitted = obj.computeQuadrature(obj.meshUnfitted.geometryType);
            
            posGP_iso_unfitted = obj.computePosGP(obj.meshUnfitted.coord_iso_per_cell,interpolation_unfitted,quadrature_unfitted);
            
            shapeValues = zeros(size(obj.meshUnfitted.connec,1),interpolation_background.nnode);
            for isubcell = 1:size(obj.meshUnfitted.connec,1) % !! VECTORIZE THIS LOOP !!
                icell = obj.meshUnfitted.cell_containing_subcell(isubcell);
                inode = obj.meshBackground.connec(icell,:);
                
                interpolation_background.computeShapeDeriv(posGP_iso_unfitted(:,:,isubcell)');
                
                djacob = obj.mapping(obj.meshUnfitted.coord(obj.meshUnfitted.connec(isubcell,:),:),interpolation_unfitted.dvolu); % !! Could be done through Geometry class?? !!
                
                F0 = (interpolation_background.shape*quadrature_unfitted.weigp')'*F1(inode)/interpolation_unfitted.dvolu;
                shapeValues(isubcell,:) = shapeValues(isubcell,:) + (interpolation_background.shape*(djacob.*quadrature_unfitted.weigp')*F0)';
            end
        end
        
        function M2 = rearrangeOutputRHS(obj,shapeValues_AllCells)
            interpolation = Interpolation.create(obj.meshBackground,'LINEAR');
            
            M2 = zeros(interpolation.npnod,1);
            for inode = 1:interpolation.nnode
                M2 = M2 + accumarray(obj.meshBackground.connec(:,inode),shapeValues_AllCells(:,inode),[interpolation.npnod,1],@sum,0);
            end
        end
        
        function itIs = isLeveSetCuttingMesh(obj)
            itIs = ~isempty(obj.meshUnfitted.backgroundCutCells);
        end
        
    end
    
    methods (Static, Access = protected)
        
        function quadrature = computeQuadrature(geometryType)
            quadrature = Quadrature.set(geometryType);
            quadrature.computeQuadrature('LINEAR');
        end
        
    end
    
    methods (Access = private)
        
        function updateMeshes(obj,unfitted)
            obj.updateBackgroundMesh(unfitted);
            obj.updateUnfittedMesh(unfitted);
        end
        
        function updateUnfittedMesh(obj,unfitted)
            obj.meshUnfitted = unfitted;
        end
        
        function updateBackgroundMesh(obj,unfitted)
            obj.meshBackground = unfitted.meshBackground;
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

