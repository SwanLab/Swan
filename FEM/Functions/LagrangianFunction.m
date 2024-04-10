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

        function fxV = evaluate(obj, xV)
            shapes = obj.interpolation.computeShapeFunctions(xV);
            nNode  = size(shapes,1);
            nGaus  = size(shapes,2);
            nF     = size(obj.fValues,2);
            nElem  = size(obj.connec,1);
            fxV = zeros(nF,1,nGaus,nElem);
            for iGaus = 1:nGaus
                for iNode = 1:nNode
                    node = (obj.connec(:,(iNode-1)*obj.ndimf+1)-1)/obj.ndimf+1;
                    Ni = shapes(iNode,iGaus);
                    fi = obj.fValues(node,:);
                    f(:,1,:) = Ni*fi';
                    fxV(:,:,iGaus,:) = squeezeParticular(fxV(:,:,iGaus,:),3) + f;
                end
            end
        end

        function fxV = sampleFunction(obj,xP,cells)
            shapes  = obj.interpolation.computeShapeFunctions(xP);
            nNode   = size(shapes,1);
            nF      = size(obj.fValues,2);
            nPoints = size(xP,2);
            fxV = zeros(nF,nPoints);
            for iF = 1:nF
                for iNode = 1:nNode
                    node = obj.mesh.connec(cells,iNode);
                    Ni = shapes(iNode,:)';
                    fi = obj.fValues(node,:);
                    f(1,:) = fi.*Ni;
                    fxV(iF,:) = fxV(iF,:) + f;
                end
            end
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
        
        function dNdx  = evaluateCartesianDerivatives(obj,xV)
            nElem = size(obj.connec,1);
            nNodeE = obj.interpolation.nnode;
            nDimE = obj.interpolation.ndime;
            nDimG = obj.mesh.ndim;
            nPoints = size(xV, 2);
            invJ  = obj.mesh.computeInverseJacobian(xV);
            deriv = obj.computeShapeDerivatives(xV);
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

        function fVR = evaluateGradient(obj, xV)
            dNdx = obj.evaluateCartesianDerivatives(xV);
            nDimf = obj.ndimf;
            nDimG = size(dNdx, 1);
            nNodeE = size(dNdx, 2);
            nPoints = size(dNdx, 3);
            nElem = size(dNdx, 4);
           
            fV = reshape(obj.fValues',[numel(obj.fValues) 1]);
            grad = zeros(nDimG,nDimf, nPoints, nElem);
            for iDimG = 1:nDimG
                for jDimf = 1:nDimf
                    for kNodeE = 1:nNodeE
                        dNdxIK = squeezeParticular(dNdx(iDimG, kNodeE,:,:),[1 2]);
                        iDofE = nDimf*(kNodeE-1)+jDimf;
                        dofs = obj.connec(:,iDofE);
                        fKJ = repmat(fV(dofs),[nPoints 1]);
                        gradIJ= dNdxIK.*fKJ;
                        grad(iDimG,jDimf,:,:) = squeezeParticular(grad(iDimG,jDimf,:,:),[1 2]) + gradIJ;
                    end
                end
            end
            fVR = reshape(grad, [nDimG*nDimf,nPoints, nElem]);
        end

        function hessianfV = sampleHessian(obj,xV)
            dNdx = obj.evaluateCartesian2Derivatives(xV);
            nDimf = obj.ndimf;
            nDimG = size(dNdx, 1);
            nNodeE = size(dNdx, 2);
            nPoints = size(dNdx, 3);
            nElem = size(dNdx, 4);
           
            fV = reshape(obj.fValues',[numel(obj.fValues) 1]);
            hessian = zeros(nDimG,nDimf, nPoints, nElem);
            for iDimG = 1:nDimG
                for jDimf = 1:nDimf
                    for kNodeE = 1:nNodeE
                        d2NdxIK = squeezeParticular(dNdx(iDimG, kNodeE,:,:),[1 2]);
                        iDofE = nDimf*(kNodeE-1)+jDimf;
                        dofs = obj.connec(:,iDofE);
                        fKJ = repmat(fV(dofs),[nPoints 1]);
                        hessianIJ= d2NdxIK.*fKJ;
                        hessian(iDimG,jDimf,:,:) = squeezeParticular(hessian(iDimG,jDimf,:,:),[1 2]) + hessianIJ;
                    end
                end
            end
            hessianfV = reshape(hessian, [nDimG*nDimf,nPoints, nElem]);
        end

        function symGrad = evaluateSymmetricGradient(obj,xV)
            grad = obj.evaluateGradient(xV);
            nDimf = obj.ndimf;
            nDims = size(grad, 1)/nDimf;
            nGaus = size(grad, 2);
            nElem = size(grad, 3);

            gradReshp = reshape(grad, [nDims,nDimf,nGaus,nElem]);
            gradT = permute(gradReshp, [2 1 3 4]);
            symGrad = 0.5*(gradReshp + gradT);
            
            symGrad =  reshape(symGrad, [nDims*nDimf,nGaus,nElem]);
        end

        function symGrad = evaluateSymmetricGradientVoigt(obj,xV)
            sgr = obj.evaluateSymmetricGradient(xV);
            symGrad = obj.obtainVoigtFormat(sgr);
        end

        function newObj = obtainVoigtFormat(obj, sgr)
            ndim = size(sgr,1);
            switch ndim
                case 4
                    newObj = obj.applyVoigt2D(sgr);
                case 9
                    newObj = obj.applyVoigt3D(sgr);
            end
        end

        function fV = applyVoigt2D(obj, fV)
            nGaus = size(fV,2);
            nElem = size(fV,3);
            fVal(1,:,:) = fV(1,:,:); % xx
            fVal(2,:,:) = fV(4,:,:); % yy
            fVal(3,:,:) = fV(2,:,:) + fV(3,:,:); % xy
            fV = reshape(fVal, [3 nGaus nElem]);
%             newObj = FGaussDiscontinuousFunction.create(fV,obj.mesh,obj.quadrature);
        end

        function fV = applyVoigt3D(obj, fV)
            nGaus = size(fV,2);
            nElem = size(fV,3);
            fVal(1,:,:) = fV(1,:,:); % xx
            fVal(2,:,:) = fV(5,:,:); % yy
            fVal(3,:,:) = fV(9,:,:); % zz
            fVal(4,:,:) = fV(2,:,:) + fV(4,:,:); % xy
            fVal(5,:,:) = fV(3,:,:) + fV(7,:,:); % xz
            fVal(6,:,:) = fV(6,:,:) + fV(8,:,:); % yz
            fV = reshape(fVal, [6 nGaus nElem]);
%             newObj = FGaussDiscontinuousFunction.create(fV,obj.mesh,obj.quadrature);
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
                        %             colorbar
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