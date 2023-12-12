classdef RHSintegrator_CutMesh < handle

    properties (Access = private)
        mesh
        quadrature
        npnod
        globalConnec
        xGauss
        fGauss

        backgroundMeshType

        xCoordsIso
        cellContainingSubcell
        subCellConnec
    end

    methods (Access = public)

        function obj = RHSintegrator_CutMesh(cParams)
            obj.init(cParams);
            obj.createQuadrature(cParams);
        end

        function rhs = compute(obj, unfFun, test)
            obj.computeSubCellConnec();
            obj.computeGaussPoints();
            obj.computeFgauss(unfFun);
            rhsCut = obj.computeElementalRHS(test);
            rhsCells = obj.assembleSubcellsInCells(rhsCut,test);
            rhs = obj.assembleIntegrand(rhsCells,test);
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

        function createQuadrature(obj,cParams)
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature(cParams.quadType);
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

        function computeFgauss(obj, unfFun)
            fG = unfFun.evaluateCutElements(obj.quadrature);
            fG = permute(fG,[2 3 1]);
            obj.fGauss = fG;
        end
        
        function rhsC = computeElementalRHS(obj,test)
            fG     = obj.fGauss;
            q      = obj.quadrature;
            dV     = obj.mesh.computeDvolume(q);
            shapes = obj.evalShapes(test);
            dofs   = size(shapes,1);
            nelem  = obj.mesh.nelem;
            int = zeros(dofs,nelem);
            for igaus = 1:obj.quadrature.ngaus
                fdv = fG(igaus,:).*dV(igaus,:);
                shape = shapes(:,:,igaus);
                int = int + bsxfun(@times,shape,fdv);
            end
            rhsC = transpose(int);
        end
        
        function rhsCells = assembleSubcellsInCells(obj,rhsCut,test)
            dofs = test.nDofsElem;
            nelem  = size(obj.globalConnec,1);
            cellNum = obj.cellContainingSubcell;
            totalInt = zeros(nelem,dofs);
            for idof = 1:dofs
                int = rhsCut(:,idof);
                intGlobal = accumarray(cellNum,int,[nelem,1],@sum,0);
                totalInt(:,idof) = totalInt(:,idof) + intGlobal;
            end
            rhsCells = totalInt;
        end

        function f = assembleIntegrand(obj,rhsCells,test)
            integrand = rhsCells;
            ndofs = test.nDofs;
            connec = test.computeDofConnectivity()';
            ndof   = size(connec,2);
            f = zeros(ndofs,1);
            for idof = 1:ndof
                int = integrand(:,idof);
                con = connec(:,idof);
                f = f + accumarray(con,int,[ndofs,1],@sum,0);
            end
        end

        function shapeAtGauss = evalShapes(obj,test)
            m.type = obj.backgroundMeshType;
            int    = Interpolation.create(m,test.order);
            int.computeShapeDeriv(obj.xGauss);
            shapeAtGauss = permute(int.shape,[1 3 2]);
        end

    end

end