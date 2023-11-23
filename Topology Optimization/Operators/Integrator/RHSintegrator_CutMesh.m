classdef RHSintegrator_CutMesh < RHSintegrator

    properties (Access = private)
        npnod
        globalConnec
        xGauss
        fGauss
        quadOrder

        backgroundMeshType

        xCoordsIso
        cellContainingSubcell
        subCellConnec
    end

    methods (Access = public)

        function obj = RHSintegrator_CutMesh(cParams)
            obj.init(cParams);
            obj.setQuadratureOrder(cParams);
            obj.createQuadrature();
        end

        function rhs = compute(obj, unfFun)
            obj.computeSubCellConnec();
            obj.computeGaussPoints();
            obj.computeFgauss(unfFun);
            rhsCut = obj.computeElementalRHS();
            rhsCells = obj.assembleSubcellsInCells(rhsCut);
            rhs = obj.assembleIntegrand(rhsCells);
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.mesh         = cParams.mesh;
            obj.npnod        = cParams.npnod;
            obj.globalConnec = cParams.globalConnec;
            obj.xCoordsIso   = cParams.xCoordsIso;
            obj.backgroundMeshType    = cParams.backgroundMeshType;
            obj.cellContainingSubcell = cParams.cellContainingSubcell;
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

        function computeFgauss(obj, unfFun)
            fG = unfFun.evaluateCutElements(obj.xGauss);
            fG = permute(fG,[2 3 1]);
            obj.fGauss = fG;
        end
        
        function rhsC = computeElementalRHS(obj)
            fG     = obj.fGauss;
            dV     = obj.computeDvolume();
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