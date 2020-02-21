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
           
            shapes = obj.computeShapes();
            dvolum = obj.computeJacobians();
            F0V = obj.interpolateFunctionInGaussPoints(F1,shapes);
            
            nelem = size(obj.meshUnfitted.connec,1);
            nnode = obj.backgroundInterp.nnode;
            shapeValues = zeros(nelem,nnode);
            for isubcell = 1:nelem % !! VECTORIZE THIS LOOP !!                
%                 shape = obj.computeShape(isubcell);
%                 weigth = obj.quadrature.weigp';
%                 djacob = obj.computeJacobian(isubcell);
%                 dvolu(1,:) = djacob*weigth;                
%                 icell  = obj.meshUnfitted.cellContainingSubcell(isubcell);
%                 inode = obj.meshBackground.connec(icell,:);
%                 Fnodes(:,1) = F1(inode);
%                 F0(1,:) = (shape)'*Fnodes;
                
                shape = shapes(:,:,isubcell);                     
                dvolu(1,:) = dvolum(isubcell,:);
                F0(1,:) = F0V(:,isubcell);

                
                shapeValues(isubcell,:) = shapeValues(isubcell,:) + (sum(dvolu.*shape.*F0,2))';
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
            posgp = zeros(size(subcell_coord,3),quadrature.ngaus,size(subcell_coord,1));
            for igaus = 1:quadrature.ngaus
                for idime = 1:size(subcell_coord,3)
                    posgp(idime,igaus,:) = subcell_coord(:,:,idime)*interpolation.shape(:,igaus);
                end
            end
        end
        
    end
    
    
    methods (Access = private)
        
      function shape = computeShape(obj,isubcell)
           xGauss = obj.unfittedQuad(:,:,isubcell);
           obj.backgroundInterp.computeShapeDeriv(xGauss);
           shape = obj.backgroundInterp.shape;
        end        
        
        function shapes = computeShapes(obj)
            xGauss = obj.unfittedQuad(:,:,:);
            obj.backgroundInterp.computeShapeDeriv(xGauss);
            shapes = obj.backgroundInterp.shape;                        
        end
        
        function computeThisQuadrature(obj)
            type = obj.meshUnfitted.geometryType;
            obj.quadrature = obj.computeQuadrature(type);
        end
        
        function computeUnfittedGaussPoints(obj)
            coord = obj.meshUnfitted.subcellIsoCoords;
            inter = obj.unfittedInterp;
            quad  = obj.quadrature;
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
            dJ = obj.mapping(coord); % !! Could be done through Geometry class?? !!
        end
        
        function dvolu = computeJacobians(obj)
           s.mesh = obj.meshUnfitted;
           g = Geometry.create(s);
           g.computeGeometry(obj.quadrature,obj.unfittedInterp)
           dvolu = g.dvolu;
        end
        
        function Fgaus = interpolateFunctionInGaussPoints(obj,F,shapes)
            connec = obj.meshUnfitted.globalConnec;
            nCell = obj.meshUnfitted.nelem;
            nnode = obj.meshBackground.nnode;         
            ngaus = size(shapes,2);
            Fgaus = zeros(ngaus,nCell);
             for inode = 1:nnode
                nodes  = connec(:,inode);
                Fnodes = F(nodes,1);
                for igaus = 1:ngaus
                   shape = squeeze(shapes(inode,igaus,:)); 
                   Fgaus(igaus,:) = Fgaus(igaus,:) + (shape.*Fnodes)';
                end
             end
        end
        
        
        function djacob = mapping(obj,coord)
            % !! PERFORM THROUGH GEOMETRY CLASS OR EXTRACT "mapping" CAPACITY FROM
            % GEOMETRY CLASS !!
            
            N_points = size(coord,1);
            switch N_points
                case 2
                    v = diff(coord);
                    L = norm(v);
                    djacob = L;
                case 3
                    if size(coord,2) == 2
                        coord = [coord, zeros(N_points,1)];
                    end
                    v1 = diff(coord([1 2],:));
                    v2 = diff(coord([1 3],:));
                    A = 0.5*norm(cross(v1,v2));
                    djacob = A;
                case 4
                    v1 = diff(coord([1 2],:));
                    v2 = diff(coord([1 3],:));
                    v3 = diff(coord([1 4],:));
                    V = (1/6)*det([v1;v2;v3]);
                    djacob = V;
            end
            isoDv  = obj.unfittedInterp.isoDv;
            djacob = djacob/isoDv;            
        end        
        
    end
    
end

