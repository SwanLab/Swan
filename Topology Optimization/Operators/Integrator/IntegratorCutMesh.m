classdef IntegratorCutMesh < Integrator
    
    properties (GetAccess = public, SetAccess = protected)
        cutMesh
        backgroundMesh
    end
    
    properties (Access = private)
        backgroundInterp
        unfittedInterp
        quadrature
        unfittedQuad
        shapes
        cutShapes
    end
    
    methods (Access = public)
        
        function obj = IntegratorCutMesh(cParams)
            obj.init(cParams);
            obj.cutMesh = obj.mesh;
            obj.backgroundMesh = obj.cutMesh.getBackgroundMesh();
        end
        
        function A = integrate(obj,F)
            obj.initShapes();
            obj.computeElementalRHS(F);
            obj.assembleSubcellsInCells();
            A = obj.assembleIntegrand();
        end
        
    end
    
    methods (Access = private)
        
        function initShapes(obj)
            nelem = obj.backgroundMesh.nelem;
            cNelem = obj.cutMesh.nelem;
            nnode = obj.backgroundMesh.nnode;
            obj.shapes = zeros(nelem,nnode);
            obj.cutShapes = zeros(cNelem,nnode);
        end
        
        function computeElementalRHS(obj,F1)
            int = obj.cutShapes;
            obj.createBackgroundInterpolation();
            obj.createUnfittedInterpolation();
            obj.computeThisQuadrature();
            obj.computeUnfittedGaussPoints();
            
            nelem = obj.cutMesh.nelem;
            nnode = obj.backgroundInterp.nnode;
            for isubcell = 1:nelem
                shape = obj.computeShape(isubcell);
                djacob = obj.computeJacobian(isubcell);
                icell  = obj.cutMesh.cellContainingSubcell(isubcell);
                inode = obj.backgroundMesh.connec(icell,:);
                weight = obj.quadrature.weigp';
                dvolu  = obj.unfittedInterp.dvolu;
                
                F0 = (shape*weight)'*F1(inode)/dvolu;
                
                int(isubcell,:) = int(isubcell,:) + (shape*(djacob.*weight)*F0)';
            end
            obj.cutShapes = int;
        end
        
        function assembleSubcellsInCells(obj)
            nnode = obj.backgroundMesh.nnode;
            nelem = obj.backgroundMesh.nelem;
            cellNum = obj.cutMesh.cellContainingSubcell;
            totalInt = obj.shapes;
            
            for iNode = 1:nnode
                int = obj.cutShapes(:,iNode);
                intGlobal  = accumarray(cellNum,int,[nelem,1],@sum,0);
                totalInt(:,iNode) = totalInt(:,iNode) + intGlobal;
            end
            obj.shapes = totalInt;
        end
        
        function f = assembleIntegrand(obj)
            integrand = obj.shapes;
            npnod  = obj.backgroundMesh.npnod;
            nnode  = obj.backgroundMesh.nnode;
            connec = obj.backgroundMesh.connec;
            f = zeros(npnod,1);
            for inode = 1:nnode
                int = integrand(:,inode);
                con = connec(:,inode);
                f = f + accumarray(con,int,[npnod,1],@sum,0);
            end
        end
        
        function itIs = isLeveSetCuttingMesh(obj)
            itIs = ~isempty(obj.cutMesh.backgroundCutCells);
        end
        
    end
    
    methods (Static, Access = private)
        
        function posgp = computePosGP(coord,interpolation,quadrature)
            interpolation.computeShapeDeriv(quadrature.posgp);
            posgp = zeros(quadrature.ngaus,size(coord,3),size(coord,1));
            for igaus = 1:quadrature.ngaus
                for idime = 1:size(coord,3)
                    posgp(igaus,idime,:) = coord(:,:,idime)*interpolation.shape(:,igaus);
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
            type = obj.cutMesh.geometryType;
            obj.quadrature = obj.computeQuadrature(type);
        end
        
        function computeUnfittedGaussPoints(obj)
            coord = obj.cutMesh.subcellIsoCoords;
            inter = obj.unfittedInterp;
            quad  =  obj.quadrature;
            quadU = obj.computePosGP(coord,inter,quad);
            obj.unfittedQuad = quadU;
        end
        
        function createBackgroundInterpolation(obj)
            mesh = obj.backgroundMesh;
            int = Interpolation.create(mesh,'LINEAR');
            obj.backgroundInterp = int;
        end
        
        function createUnfittedInterpolation(obj)
            mesh = obj.cutMesh;
            int = Interpolation.create(mesh,'LINEAR');
            obj.unfittedInterp = int;
        end
        
        function dJ = computeJacobian(obj,isubcell)
            connec = obj.cutMesh.connec(isubcell,:);
            coord  = obj.cutMesh.coord(connec,:);
            dvolu = obj.unfittedInterp.dvolu;
            dJ = obj.mapping(coord,dvolu); % !! Could be done through Geometry class?? !!
        end
        
    end
    
end

