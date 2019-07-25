classdef IntegratorUnfitted < Integrator
    
    properties (GetAccess = public, SetAccess = protected)
        meshUnfitted
        meshBackground
    end
    
    properties (Access = private)
        backgroundInterp        
        unfittedInterp
        quadrature
        unfittedQuad        
    end    
    
    methods (Access = public, Abstract)
        integrate(obj)
    end
    
    methods (Access = public)
        
        function obj = IntegratorUnfitted(cParams)
            obj.init(cParams);
            obj.meshUnfitted   = obj.mesh;            
            obj.meshBackground = obj.mesh.meshBackground;                   
        end
        
        function A = integrateUnfittedMesh(obj,F,meshUnfitted)
            if exist('meshUnfitted','var')
                %obj.updateMeshes(meshUnfitted);
             %   obj.updateBackgroundMesh();
             %   obj.updateUnfittedMesh();   
            end
            A = obj.integrate(F);
        end
        

        
    end
    
    methods (Static, Access = public)

        
    end
    
    methods (Access = protected)
        
        function shapeValues = evaluateCutShapes(obj,F1)
            obj.createBackgroundInterpolation();
            obj.createUnfittedInterpolation();
            obj.computeThisQuadrature();
            obj.computeUnfittedGaussPoints();
            
            
            nelem = size(obj.meshUnfitted.connec,1);
            nnode = obj.backgroundInterp.nnode;
            shapeValues = zeros(nelem,nnode);
            for isubcell = 1:nelem % !! VECTORIZE THIS LOOP !!
                
                shape = obj.computeShape(isubcell);
                
                djacob = obj.computeJacobian(isubcell);
                icell  = obj.meshUnfitted.cellContainingSubcell(isubcell);

                inode = obj.meshBackground.connec(icell,:);
                
                weigth = obj.quadrature.weigp';
                dvolu  = obj.unfittedInterp.dvolu;
                
                F0 = (shape*weigth)'*F1(inode)/dvolu;
                
                shapeValues(isubcell,:) = shapeValues(isubcell,:) + (shape*(djacob.*weigth)*F0)';
            end
        end      
        

        
        function M2 = rearrangeOutputRHS(obj,shapes)
            npnod = obj.meshBackground.npnod;
            nnode = obj.meshBackground.nnode;
            
            M2 = zeros(npnod,1);
            for inode = 1:nnode
                M2 = M2 + accumarray(obj.meshBackground.connec(:,inode),shapes(:,inode),[npnod,1],@sum,0);
            end
        end
        
        function itIs = isLeveSetCuttingMesh(obj)
            itIs = ~isempty(obj.meshUnfitted.backgroundCutCells);
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
    
    
    methods (Access = private)
        
        function shape = computeShape(obj,isubcell)
           xGauss = obj.unfittedQuad(:,:,isubcell)';
           obj.backgroundInterp.computeShapeDeriv(xGauss);
           shape = obj.backgroundInterp.shape;
        end
        
        function computeThisQuadrature(obj)
            type = obj.meshUnfitted.geometryType;
            obj.quadrature = obj.computeQuadrature(type);
        end
        
        function computeUnfittedGaussPoints(obj)
           coord = obj.meshUnfitted.subcellIsoCoords;
           inter = obj.unfittedInterp;
           quad  =  obj.quadrature;
           quadU = obj.computePosGP(coord,inter,quad);
           obj.unfittedQuad = quadU;
        end
        
        function createBackgroundInterpolation(obj)
            mesh = obj.meshBackground;
            int = Interpolation.create(mesh,'LINEAR');            
            obj.backgroundInterp = int;
        end
        
       function createUnfittedInterpolation(obj)
           mesh = obj.meshUnfitted;
           int = Interpolation.create(mesh,'LINEAR');
           obj.unfittedInterp = int;            
        end
        
        function dJ = computeJacobian(obj,isubcell)
            connec = obj.meshUnfitted.connec(isubcell,:);
            coord  = obj.meshUnfitted.coord(connec,:);
            dvolu = obj.unfittedInterp.dvolu;
            dJ = obj.mapping(coord,dvolu); % !! Could be done through Geometry class?? !!            
        end
        
    end
    
end

