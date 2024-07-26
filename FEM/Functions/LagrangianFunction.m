classdef LagrangianFunction < FeFunction

    properties (GetAccess = public, SetAccess = private)
        nDofs
        nDofsElem
        connec
    end

    properties (Access = private)
        interpolation
        coord
    end

    methods (Access = public)

        function obj = LagrangianFunction(cParams)
            obj.init(cParams);
            obj.createInterpolation();
            obj.createDOFCoordConnec();
        end

        %%% Function %%%
        function fxV = evaluate(obj, xV)
            shapes  = obj.interpolation.computeShapeFunctions(xV);
            nNodesE  = size(shapes,1);
            nPoints = size(shapes,2);
            nDimf      = size(obj.fValues,2);
            nElem   = size(obj.connec,1);
            fxV = zeros(nDimf,nPoints,nElem);
            for iDimf = 1:nDimf
                for iNodeE = 1:nNodesE
                    dofE = (iNodeE-1)*nDimf + iDimf;
                    dofs = (obj.connec(:,dofE));
                    nodes = (dofs - iDimf)/nDimf + 1;
                    Ni = shapes(iNodeE,:);
                    fi = obj.fValues(nodes,iDimf);
                    f(1,:,:) = (fi*Ni)';
                    fxV(iDimf,:,:) = fxV(iDimf,:,:)  + f;
                end
            end
        end

        function fxP = sampleFunction(obj,xP,cells)
            shapes  = obj.interpolation.computeShapeFunctions(xP);
            nNodesE   = size(shapes,1);
            nPoints = size(xP,2);
            nDimf      = size(obj.fValues,2);
            fxP = zeros(nDimf,1,nPoints);
            for iDimf = 1:nDimf
                for iNodeE = 1:nNodesE
                    dofE = (iNodeE-1)*nDimf + iDimf;
                    dofs = (obj.connec(cells,dofE));
                    nodes = (dofs - iDimf)/nDimf + 1;
                    Ni = shapes(iNodeE,:);
                    fi = obj.fValues(nodes,iDimf);
                    f(1,1,:) = Ni.*fi';
                    fxP(iDimf,:,:) = fxP(iDimf,:,:)  + f;
                end
            end
        end

        %%% Jacobian %%%
        function dNdx  = evaluateCartesianDerivatives(obj,xV)
            nElem = size(obj.connec,1);
            nNodeE = obj.interpolation.nnode;
            nDimE = obj.interpolation.ndime;
            nDimG = obj.mesh.ndim;
            nPoints = size(xV, 2);
            invJ  = obj.mesh.computeInverseJacobian(xV);
            deriv = obj.interpolation.computeShapeDerivatives(xV);
            dShapes  = zeros(nDimG,nNodeE,nPoints,nElem);
            for iDimG = 1:nDimG
                for kNodeE = 1:nNodeE
                    for jDimE = 1:nDimE
                        invJ_IJ   = invJ(iDimG,jDimE,:,:);
                        dShapes_JK = deriv(jDimE,kNodeE,:);
                        dShapes_KI   = pagemtimes(invJ_IJ,dShapes_JK);
                        dShapes(iDimG,kNodeE,:,:) = dShapes(iDimG,kNodeE,:,:) + dShapes_KI;
                    end
                end
            end
            dNdx = dShapes;
        end

        function dNdx  = sampleCartesianDerivatives(obj,xP,cells)
            nNodeE = obj.interpolation.nnode;
            nDimE = obj.interpolation.ndime;
            nDimG = obj.mesh.ndim;
            nPoints = size(xP, 2);
            invJ  = obj.mesh.sampleInverseJacobian(xP,cells);
            deriv = obj.interpolation.computeShapeDerivatives(xP);
            dShapes  = zeros(nDimG,nNodeE,nPoints);
            for iDimG = 1:nDimG
                for kNodeE = 1:nNodeE
                    for jDimE = 1:nDimE
                        invJ_IJ   = invJ(iDimG,jDimE,:);
                        dShapes_JK = deriv(jDimE,kNodeE,:);
                        dShapes_KI   = pagemtimes(invJ_IJ,dShapes_JK);
                        dShapes(iDimG,kNodeE,:) = dShapes(iDimG,kNodeE,:) + dShapes_KI;
                    end
                end
            end
            dNdx = dShapes;
        end 

        function dfxV = evaluateGradient(obj, xV)
            dNdx = obj.evaluateCartesianDerivatives(xV);
            nDimf = obj.ndimf;
            nDimG = size(dNdx, 1);
            nNodeE = size(dNdx, 2);
            nPoints = size(dNdx, 3);
            nElem = size(dNdx, 4);
           
            fV = reshape(obj.fValues',[numel(obj.fValues) 1]);
            grad = zeros(nDimG, nDimf, nPoints, nElem);
            for iDimG = 1:nDimG
                for jF = 1:nDimf
                    for kNodeE = 1:nNodeE
                        dNdxIK = squeezeParticular(dNdx(iDimG, kNodeE,:,:),[1 2]);
                        dofE = (kNodeE-1)*nDimf + jF;
                        dofs = obj.connec(:,dofE);
                        fKJ = repmat(fV(dofs),[1 nPoints]);
                        gradIJ= dNdxIK.*fKJ';
                        grad(iDimG,jF,:,:) = squeezeParticular(grad(iDimG,jF,:,:),[1 2]) + gradIJ;
                    end
                end
            end
            dfxV = reshape(grad, [nDimG*nDimf,nPoints, nElem]);
        end

        function dfxP = sampleGradient(obj,xP,cells)
            dNdx = obj.sampleCartesianDerivatives(xP,cells);
            nDimf = obj.ndimf;
            nDimG = size(dNdx, 1);
            nNodeE = size(dNdx, 2);
            nPoints = size(dNdx, 3);
           
            fV = reshape(obj.fValues',[numel(obj.fValues) 1]);
            grad = zeros(nDimG,nDimf,nPoints);
            for iDimG = 1:nDimG
                for jDimf = 1:nDimf
                    for kNodeE = 1:nNodeE
                        dNdxIK = squeezeParticular(dNdx(iDimG, kNodeE,:),[1 2]);
                        iDofE = nDimf*(kNodeE-1)+jDimf;
                        dofs = obj.connec(cells,iDofE);
                        fKJ = fV(dofs);
                        gradIJ = dNdxIK.*fKJ;
                        grad(iDimG,jDimf,:) = squeezeParticular(grad(iDimG,jDimf,:),[1 2]) + gradIJ;
                    end
                end
            end
            dfxP = reshape(grad, [nDimG*nDimf,nPoints]);
        end

       %%% Hessian %%%
        function d2Ndx = evaluateCartesianSecondDerivatives(obj,xV)
            nElem = size(obj.connec,1);
            nNodeE = obj.interpolation.nnode;
            nDimE = obj.interpolation.ndime;
            nDimG = obj.mesh.ndim;
            nPoints = size(xV, 2);
            invJ  = obj.mesh.computeInverseJacobian(xV);
            deriv2 = obj.interpolation.computeShapeSecondDerivatives(xV);
            ddShapes  = zeros(nDimG,nNodeE,nPoints,nElem);
            for iDimG = 1:nDimG
                for kNodeE = 1:nNodeE
                    for jDimE = 1:nDimE
                        invJ_IJ   = invJ(iDimG,jDimE,:,:);
                        ddShapes_JK = deriv2(jDimE,kNodeE,:);
                        ddShapes_KI   = pagemtimes(invJ_IJ,ddShapes_JK);
                        ddShapes(iDimG,kNodeE,:,:) = ddShapes(iDimG,kNodeE,:,:) + ddShapes_KI;
                    end
                end
            end
            d2Ndx = ddShapes;
        end

        function d2Ndx = sampleCartesianSecondDerivatives(obj,xP,cells)
            nNodeE = obj.interpolation.nnode;
            nDimE = obj.interpolation.ndime;
            nDimG = obj.mesh.ndim;
            nPoints = size(xP, 2);
            invJ  = obj.mesh.sampleInverseJacobian(xP,cells);
            deriv2 = obj.interpolation.computeShapeSecondDerivatives(xP);
            ddShapes  = zeros(nDimG,nNodeE,nPoints);
            for iDimG = 1:nDimG
                for kNodeE = 1:nNodeE
                    for jDimE = 1:nDimE
                        invJ_IJ   = invJ(iDimG,jDimE,:);
                        ddShapes_JK = deriv2(jDimE,kNodeE,:);
                        ddShapes_KI   = pagemtimes(invJ_IJ,ddShapes_JK);
                        ddShapes(iDimG,kNodeE,:) = ddShapes(iDimG,kNodeE,:) + ddShapes_KI;
                    end
                end
            end
            d2Ndx = ddShapes;
        end

        function d2fxV = evaluateHessian(obj,xV)
            d2Ndx = obj.evaluateCartesianSecondDerivatives(xV);
            nDimf = obj.ndimf;
            nDimG = size(d2Ndx, 1);
            nNodeE = size(d2Ndx, 2);
            nPoints = size(d2Ndx, 3);
            nElem = size(d2Ndx, 4);
           
            fV = reshape(obj.fValues',[numel(obj.fValues) 1]);
            hessian = zeros(nDimG, nDimf, nPoints, nElem);
            for iDimG = 1:nDimG
                for jDimf = 1:nDimf
                    for kNodeE = 1:nNodeE
                        d2NdxIK = squeezeParticular(d2Ndx(iDimG, kNodeE,:,:),[1 2]);
                        dofE = (kNodeE-1)*nDimf + jDimf;
                        dofs = obj.connec(:,dofE);
                        fKJ = repmat(fV(dofs),[1 nPoints]);
                        hessianIJ= d2NdxIK.*fKJ';
                        hessian(iDimG,jDimf,:,:) = squeezeParticular(hessian(iDimG,jDimf,:,:),[1 2]) + hessianIJ;
                    end
                end
            end
            d2fxV = reshape(hessian, [nDimG*nDimf,nPoints, nElem]);
        end

        function d2fxP = sampleHessian(obj,xP,cells)
            d2Ndx = obj.sampleCartesianSecondDerivatives(xP,cells);
            nDimf = obj.ndimf;
            nDimG = size(d2Ndx, 1);
            nNodeE = size(d2Ndx, 2);
            nPoints = size(d2Ndx, 3);

            fV = reshape(obj.fValues',[numel(obj.fValues) 1]);
            hessian = zeros(nDimG,nDimf,nPoints);
            for iDimG = 1:nDimG
                for jDimf = 1:nDimf
                    for kNodeE = 1:nNodeE
                        d2NdxIK = squeezeParticular(d2Ndx(iDimG, kNodeE,:),[1 2]);
                        iDofE = nDimf*(kNodeE-1)+jDimf;
                        dofs = obj.connec(cells,iDofE);
                        fKJ = fV(dofs);
                        hessianIJ= d2NdxIK.*fKJ;
                        hessian(iDimG,jDimf,:) = squeezeParticular(hessian(iDimG,jDimf,:),[1 2]) + hessianIJ;
                    end
                end
            end
            d2fxP = reshape(hessian, [nDimG*nDimf,nPoints]);
        end
        
        function c = getCoord(obj)
            c = obj.coord;
        end
        
        function c = getConnec(obj)
            c = obj.connec;
        end

        function N = computeShapeFunctions(obj, xV)
            N = obj.interpolation.computeShapeFunctions(xV);
        end

        function dN = computeShapeDerivatives(obj, xV)
            dN = obj.interpolation.computeShapeDerivatives(xV);
        end

        function ddN = computeShapeSecondDerivatives(obj,xV)
            ddN = obj.interpolation.computeShapeSecondDerivatives(xV);
        end
        
        function ord = orderTextual(obj)
            switch obj.order
                case 'P0'
                    ord = 'CONSTANT';
                case 'P1'
                    ord = 'LINEAR';
                case 'P2'
                    ord = 'QUADRATIC';
                case 'P3'
                    ord = 'CUBIC';
            end
        end

        function ord = getOrderNum(obj)
            ord = str2double(obj.order(end));
        end

        function plot(obj) % 2D domains only
            if  strcmp(obj.order,'LINEAR')
                switch obj.mesh.type
                case {'TRIANGLE','QUAD'}
                    x = obj.coord(:,1);
                    y = obj.coord(:,2);
                    figure(1000)
                    for idim = 1:obj.ndimf
                        subplot(1,obj.ndimf,idim);
                        z = obj.fValues(:,idim);
                        a = trisurf(obj.connec,x,y,z);
                        view(0,90)
                        % colorbar
                        shading interp
                        a.EdgeColor = [0 0 0];
                        title(['dim = ', num2str(idim)]);
                    end
                case 'LINE'
                    x = obj.mesh.coord(:,1);
                    y = obj.fValues;
                    figure(1000)
                    plot(x,y)
                end
            else
                pl = LagrangianPlotter();
                s.func = obj;
                s.mesh = obj.mesh;
                s.interpolation = obj.interpolation;
                pl.plot(s);
            end
        end

        function dofConnec = computeDofConnectivity(obj)
            conne  = obj.connec;
            nDimf  = obj.ndimf;
            nNode  = size(conne, 2);
            nDofsE = nNode*nDimf;
            dofsElem  = zeros(nDofsE,size(conne,1));
            for iNode = 1:nNode
                for iUnkn = 1:nDimf
                    idofElem   = nDimf*(iNode - 1) + iUnkn;
                    globalNode = conne(:,iNode);
                    idofGlobal = nDimf*(globalNode - 1) + iUnkn;
                    dofsElem(idofElem,:) = idofGlobal;
                end
            end
            dofConnec = dofsElem;
        end

        function dof = getDofsFromCondition(obj, condition)
            nodes = condition(obj.coord);
            iNode = find(nodes==1);
            dofElem = repmat(1:obj.ndimf, [length(iNode) 1]);
            dofMat = obj.ndimf*(iNode - 1) + dofElem;
            dof = sort(dofMat(:));
        end

        function print(obj, filename, software)
            if nargin == 2; software = 'Paraview'; end
%             sF.fValues = obj.fValues;
%             sF.mesh = obj.mesh;
%             p1 = P1Function(sF);
            s.mesh = obj.mesh;
            s.fun = {obj};
            s.type = software;
            s.filename = filename;
            p = FunctionPrinter.create(s);
            p.print();
        end

        function [res, pformat] = getDataToPrint(obj)
            switch obj.order
                case 'P0'
                    q = Quadrature.set(obj.mesh.type);
                    q.computeQuadrature('LINEAR');
                    nElem = size(obj.mesh.connec, 1);
                    nGaus = q.ngaus;
        
                    s.nDimf   = obj.ndimf;
                    s.nData   = nElem*nGaus;
                    s.nGroup  = nElem;
                    s.fValues = obj.getFormattedP0FValues();
                    fps = FunctionPrintingSettings(s);
                    [res, pformat] = fps.getDataToPrint();

                otherwise
                    nNods = size(obj.fValues, 1);
                    s.nDimf   = obj.ndimf;
                    s.nData   = nNods;
                    s.nGroup  = nNods;
                    s.fValues = obj.fValues;
                    fps = FunctionPrintingSettings(s);
                    [res, pformat] = fps.getDataToPrint();
            end
        end
        
        function v = computeL2norm(obj)
            s.type     = 'ScalarProduct';
            s.quadType = 'QUADRATIC';
            s.mesh     = obj.mesh;
            int = Integrator.create(s);
            ff  = int.compute(obj,obj);
            v   = sqrt(ff);
        end

        function fdivF = computeFieldTimesDivergence(obj,xV)
            fG  = obj.evaluate(xV);
            dfG = obj.computeDivergence(xV);
            fdivFG = bsxfun(@times,dfG.fValues,fG);
            s.quadrature = xV;
            s.mesh       = obj.mesh;
            s.fValues    = fdivFG;
            fdivF = FGaussDiscontinuousFunction(s);
        end

        function divF = computeDivergence(obj,xV)
            dNdx = obj.evaluateCartesianDerivatives(xV);
            fV = obj.fValues;
            nodes = obj.mesh.connec;
            nNode = obj.mesh.nnodeElem;
            nDim  = obj.mesh.ndim;
            nGaus = size(xV,2);
            divV = zeros(nGaus,obj.mesh.nelem);
            for igaus = 1:nGaus
                for kNode = 1:nNode
                    nodeK = nodes(:,kNode);
                    for rDim = 1:nDim
                        dNkr = squeeze(dNdx(rDim,kNode,igaus,:));
                        fkr = fV(nodeK,rDim);
                        int(1,:) = dNkr.*fkr;
                        divV(igaus,:) = divV(igaus,:) + int;
                    end
                end
            end
            s.quadrature = xV;
            s.mesh       = obj.mesh;
            s.fValues(1,:,:) = divV;
            divF = FGaussDiscontinuousFunction(s);
        end

        function fFine = refine(obj,m,mFine)
            fNodes  = obj.fValues;
            fEdges  = obj.computeFunctionInEdges(m, fNodes);
            fAll    = [fNodes;fEdges];
            s.mesh    = mFine;
            s.fValues = fAll;
            s.order   = 'P1';
            fFine = LagrangianFunction(s);
        end

        function f = copy(obj)
            f = obj.create(obj.mesh,obj.ndimf,obj.order);
            f.fValues = obj.fValues;
        end

        function f = normalize(obj,type,epsilon)
            switch type
                case 'L2'
                    fNorm = Norm.computeL2(obj.mesh,obj);
                case 'H1'
                    fNorm = Norm.computeH1(obj.mesh,obj,epsilon);
            end
            f = obj.create(obj.mesh,obj.ndimf,obj.order);
            f.fValues = obj.fValues/sqrt(fNorm);
        end

        % Operator overload

        function s = plus(obj1,obj2)
            if isa(obj1, 'LagrangianFunction')
                res = copy(obj1);
                val1 = obj1.fValues;
            else
                val1 = obj1;
            end
            if isa(obj2, 'LagrangianFunction')
                res = copy(obj2);
                val2 = obj2.fValues;
            else
                val2 = obj2;
            end

            res.fValues = val1 + val2;
            s = res;
        end

        function s = minus(obj1,obj2)
            if isa(obj1, 'LagrangianFunction')
                res = copy(obj1);
                val1 = obj1.fValues;
            else
                val1 = obj1;
            end
            if isa(obj2, 'LagrangianFunction')
                res = copy(obj2);
                val2 = obj2.fValues;
            else
                val2 = obj2;
            end

            res.fValues = val1 - val2;
            s = res;
        end

        function r = uminus(a)
            r = copy(a);
            r.fValues = -a.fValues;
        end

        function s = times(obj1,obj2)
            res = copy(obj1);
            res.fValues = obj1.fValues .* obj2.fValues;
            s = res;
        end

        function s = power(f,b)
            res = copy(f);
            res.fValues = f.fValues .^ b;
            s = res;
        end

        function s = rdivide(f,b)
            res = copy(f);
            res.fValues = f.fValues ./ b;
            s = res;
        end

        function s = mrdivide(f,b)
            res = copy(f);
            res.fValues = f.fValues ./ b;
            s = res;
        end
        
    end

    methods (Access = public, Static)

        function pL = create(mesh, ndimf, ord)
            s.mesh    = mesh;
            s.order   = ord;
            s.ndimf   = ndimf;
            c = DOFsComputer(s);
            c.computeDofs();
            s.fValues = zeros(c.getNumberDofs()/ndimf,ndimf);
            pL = LagrangianFunction(s);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh    = cParams.mesh;
            obj.fValues = cParams.fValues;
            obj.ndimf   = size(cParams.fValues,2);
            obj.order   = cParams.order;
        end

        function createInterpolation(obj)
            type = obj.mesh.type;
            obj.interpolation = Interpolation.create(type,obj.orderTextual());
            obj.nDofsElem = obj.ndimf*obj.interpolation.nnode;
        end

        function createDOFCoordConnec(obj)
            s.mesh          = obj.mesh;
            s.interpolation = obj.interpolation;
            s.order         = obj.order;
            s.ndimf         = obj.ndimf;
            c = DOFsComputer(s);
            c.computeDofs();
            c.computeCoord();
            obj.coord  = c.getCoord();
            obj.connec = c.getDofs();
            obj.nDofs = c.getNumberDofs();
        end

        function f = computeFunctionInEdges(obj,m,fNodes)
            s.edgeMesh = m.computeEdgeMesh();
            s.fNodes   = fNodes;
            eF         = EdgeFunctionInterpolator(s);
            f = eF.compute();
        end
        function fM = getFormattedP0FValues(obj)
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature('LINEAR');
            fV = obj.evaluate(q.posgp);
            nGaus   = q.ngaus;
            nComp   = obj.ndimf;
            nElem   = size(obj.mesh.connec, 1);
            fM  = zeros(nGaus*nElem,nComp);
            for iStre = 1:nComp
                for iGaus = 1:nGaus
                    rows = linspace(iGaus,(nElem - 1)*nGaus + iGaus,nElem);
                    fM(rows,iStre) = fV(iStre,iGaus,:);
                end
            end
        end

    end

end