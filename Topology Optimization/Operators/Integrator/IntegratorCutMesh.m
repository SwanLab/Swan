classdef IntegratorCutMesh < Integrator
    
    properties (Access = private)
        Fnodal
        Fprojection
        quadrature
        shapes
        RHScells
        RHScellsCut
        backgroundMesh
    end
    
    methods (Access = public)
        
        function obj = IntegratorCutMesh(cParams)
            obj.init(cParams);
            obj.backgroundMesh = cParams.meshBackground;
        end
        
        function A = integrate(obj,Fnodal,quad)
            obj.Fnodal = Fnodal;
            if nargin == 3
                obj.quadrature = quad;
            else
                type = obj.mesh.geometryType;
                obj.quadrature = obj.computeQuadrature(type);
            end
            obj.computeElementalRHS();
            obj.assembleSubcellsInCells();
            A = obj.assembleIntegrand();
        end
        
    end
    
    methods (Access = private)
        
        function computeElementalRHS(obj)
            obj.computeShapesInUnfittedGaussPoints();
            int = obj.integrateFwithShapeFunction();
            obj.RHScellsCut = int;
        end
        
        function computeShapesInUnfittedGaussPoints(obj)
            xG    = obj.computeUnfittedGaussPoints();
            shape = obj.createShapes(obj.backgroundMesh,xG);
            obj.shapes = shape;
        end
        
        function xGauss = computeUnfittedGaussPoints(obj)
            quad  = obj.quadrature;
            coord = obj.mesh.subcellIsoCoords;
            coord = permute(coord,[1 3 2]);
            shape = obj.createShapes(obj.mesh,quad.posgp);
            nDime = size(coord,2);
            nNode = obj.mesh.nnode;
            nElem = obj.mesh.nelem;
            nGaus = quad.ngaus;
            xGauss = zeros(nGaus,nElem,nDime);
            for kNode = 1:nNode
                shapeKJ(:,1) = shape(kNode,:);
                xKJ(1,:,:) = coord(:,:,kNode);
                xG = bsxfun(@times,shapeKJ,xKJ);
                xGauss = xGauss + xG;
            end
            xGauss = permute(xGauss,[3 1 2]);
        end
        
        function Fproj = integrateFwithShapeFunction(obj)
            Fgauss = obj.interpolateFunctionInGaussPoints();
            dvolume = obj.mesh.computeDvolume(obj.quadrature);
            fdV = (Fgauss.*dvolume);
            nelem = obj.mesh.nelem;
            nnode = obj.backgroundMesh.nnode;
            Fproj = zeros(nnode,nelem);
            for igaus = 1:obj.quadrature.ngaus
                fdv = fdV(igaus,:);
                shape = obj.shapes(:,:,igaus);
                Fproj = Fproj + bsxfun(@times,shape,fdv);
            end
            Fproj = Fproj';
        end
        
        function Fgaus = interpolateFunctionInGaussPoints(obj)
            nCell  = obj.mesh.nelem;
            ngaus = obj.quadrature.ngaus;
            Fgaus = zeros(ngaus,nCell);
            Fnodes = obj.computeFinNodesPerElement();
            for igaus = 1:ngaus
                shape = obj.shapes(:,:,igaus);
                int = sum(shape.*Fnodes,1);
                Fgaus(igaus,:) = Fgaus(igaus,:) + int;
            end
        end
        
        function Fnodes = computeFinNodesPerElement(obj)
            connec = obj.mesh.globalConnec;
            nCell  = obj.mesh.nelem;
            nnode  = obj.backgroundMesh.nnode;
            Fnodes = zeros(nnode,nCell);
            for inode = 1:nnode
                nodes  = connec(:,inode);
                Fnodes(inode,:) = obj.Fnodal(nodes,1);
            end
        end
        
        function assembleSubcellsInCells(obj)
            nnode = obj.backgroundMesh.nnode;
            nelem = obj.backgroundMesh.nelem;
            cellNum = obj.mesh.cellContainingSubcell;
            totalInt = zeros(nelem,nnode);
            for iNode = 1:nnode
                int = obj.RHScellsCut(:,iNode);
                intGlobal  = accumarray(cellNum,int,[nelem,1],@sum,0);
                totalInt(:,iNode) = totalInt(:,iNode) + intGlobal;
            end
            obj.RHScells = totalInt;
        end
        
        function f = assembleIntegrand(obj)
            integrand = obj.RHScells;
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
        
    end
    
    methods (Access = private, Static)
        
        function shapes = createShapes(mesh,xG)
            int = Interpolation.create(mesh,'LINEAR');
            int.computeShapeDeriv(xG);
            shapes = permute(int.shape,[1 3 2]);
        end
        
    end
    
end

