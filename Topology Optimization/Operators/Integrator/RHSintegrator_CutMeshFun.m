classdef RHSintegrator_CutMeshFun < handle

    properties (Access = private)
        npnod
        mesh
        globalConnec
        xGauss
        fGauss
        quadOrder
        quadrature

        backgroundMeshType

        xCoordsIso
        cellContainingSubcell
        subCellConnec

        bgMesh % delete
    end

    methods (Access = public)

        % Via Integrator_Simple + Integrator
        function obj = RHSintegrator_CutMeshFun(cParams)
            obj.init(cParams);
            obj.createQuadrature();
        end

        function rhs = compute(obj, fNodal)
            obj.computeSubCellConnec();
            obj.computeGaussPoints();
            obj.computeFgauss(fNodal);
            rhsCut = obj.computeElementalRHS();
            rhsCells = obj.assembleSubcellsInCells(rhsCut);
            rhs = obj.assembleIntegrand(rhsCells);
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.mesh         = cParams.mesh;
            obj.npnod        = cParams.npnod;
            obj.quadOrder    = 'LINEAR';
            obj.globalConnec = cParams.globalConnec;

            obj.backgroundMeshType  = cParams.backgroundMeshType;

            obj.xCoordsIso   = cParams.xCoordsIso;
            obj.cellContainingSubcell = cParams.cellContainingSubcell;

            %
            obj.bgMesh = cParams.backgroundMesh;
            fV = zeros(length(obj.cellContainingSubcell),1);
            fV(obj.cellContainingSubcell) = 1;
            s.fValues = fV;
            s.mesh = obj.bgMesh;
            p0cells = P0Function(s);
    
        end

        function createQuadrature(obj)
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature('LINEAR');
            obj.quadrature = q;
        end
        
        function computeSubCellConnec(obj)
            cells = obj.cellContainingSubcell;
            obj.subCellConnec = obj.globalConnec(cells,:);
        end
        
        function computeGaussPoints(obj)
            q = obj.quadrature;
            msh.connec = obj.computeSubCellsLocalConnec();
            msh.type    = obj.mesh.type;
            s.fValues = obj.computeSubCellsLocalCoord();
            s.mesh    = msh;
            x = P1Function(s);
            obj.xGauss = x.evaluate(q.posgp);
        end
        
        function c = computeSubCellsLocalCoord(obj)
            coord = obj.xCoordsIso; 
            nDim  = size(coord,1);
            c = reshape(coord,nDim,[])';
        end
        
        function lConnec = computeSubCellsLocalConnec(obj)
            coord = obj.xCoordsIso;
            nElem = size(coord,3);
            nNode = size(coord,2);
            lConnec = reshape(1:nElem*nNode,nNode,nElem)';
        end

        function computeFgauss(obj, fNodal)
            s.fValues = fNodal;
            mmm.connec = obj.subCellConnec;
            mmm.type   = obj.backgroundMeshType;
            s.mesh = mmm;
%             s.connec = obj.subCellConnec;
%             s.type   = obj.backgroundMeshType;
%             s.mesh   = obj.mesh; % !!!
            f = P1Function(s);
            fG = f.evaluate(obj.xGauss);
            fG = permute(fG,[2 3 1]);
            obj.fGauss = fG;
        end
        
        function rhsC = computeElementalRHS(obj) % integrate@RHSintegrator
            fG     = obj.fGauss;
            dV     = obj.computeDvolume();
%             fdV    = (fG.*dV);
            shapes = obj.computeShapeFunctions();
            nnode  = size(shapes,1);
            nelem  = size(shapes,2);
            int = zeros(nnode,nelem);
            for igaus = 1:obj.quadrature.ngaus
                fdv = fG(igaus,:).*dV(igaus,:);
                shape = shapes(:, :, igaus);
                int = int + bsxfun(@times,shape,fdv);
            end
            rhsC = transpose(int);
        end
        
        function rhsCells = assembleSubcellsInCells(obj,rhsCut)
            nnode = size(obj.globalConnec,2);
            nelem = size(obj.globalConnec,1);
            cellNum = obj.cellContainingSubcell;
            totalInt = zeros(nelem,nnode);
            for iNode = 1:nnode
                int = rhsCut(:,iNode);
                intGlobal = accumarray(cellNum,int,[nelem,1],@sum,0);
                totalInt(:,iNode) = totalInt(:,iNode) + intGlobal;
            end
            rhsCells = totalInt;
        end

        function f = assembleIntegrand(obj,rhsCells)
            integrand = rhsCells;
            ndofs = obj.npnod;
            connec = obj.globalConnec;
            nnode  = size(connec,2);
            f = zeros(ndofs,1);
            for inode = 1:nnode
                int = integrand(:,inode);
                con = connec(:,inode);
                f = f + accumarray(con,int,[ndofs,1],@sum,0);
            end
        end
        
        function dV = computeDvolume(obj)
            q = obj.quadrature;
            dV = obj.mesh.computeDvolume(q);
        end
        
        function shapes = computeShapeFunctions(obj)
            m.type = obj.backgroundMeshType;
            int = Interpolation.create(m,'LINEAR');
            int.computeShapeDeriv(obj.xGauss);
            shapes = permute(int.shape,[1 3 2]);
        end

    end

end