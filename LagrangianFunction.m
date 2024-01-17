classdef LagrangianFunction < FeFunction

    properties (Access = public)
    end

    properties (Access = private)
        interpolation
        coord
        connec
    end

    methods (Access = public)

        function obj = LagrangianFunction(cParams)
            obj.init(cParams);
            obj.createInterpolation();
            obj.createDOFCoordConnec();
        end

        function fxV = evaluate(obj, xV)
            obj.interpolation.computeShapeDeriv(xV);
            shapes = obj.interpolation.shape;
            nNode  = size(shapes,1);
            nGaus  = size(shapes,2);
            nF     = size(obj.fValues,2);
            nElem  = size(obj.connec,1);
            fxV = zeros(nF,nGaus,nElem);
            for iGaus = 1:nGaus
                for iNode = 1:nNode
                    node = obj.connec(:,iNode);
                    Ni = shapes(iNode,iGaus);
                    fi = obj.fValues(node,:);
                    f(:,1,:) = Ni*fi';
                    fxV(:,iGaus,:) = fxV(:,iGaus,:) + f;
                end
            end

        end
        
        function c = getCoord(obj)
            c = obj.coord;
        end
        
        function c = getConnec(obj)
            c = obj.connec;
        end

        function N = computeShapeFunctions(obj, quad)
%             obj.mesh.computeInverseJacobian(quad,obj.interpolation);
            xV = quad.posgp;
            obj.interpolation.computeShapeDeriv(xV);
            N = obj.interpolation.shape;
        end
        
        function dNdx  = computeCartesianDerivatives(obj,quad)
            nElem = size(obj.connec,1);
            nNode = obj.interpolation.nnode;
            nDime = obj.interpolation.ndime;
            nGaus = quad.ngaus;
            invJ  = obj.mesh.computeInverseJacobian(quad,obj.interpolation);
            dShapeDx  = zeros(nDime,nNode,nElem,nGaus);
            for igaus = 1:nGaus
                dShapes = obj.interpolation.deriv(:,:,igaus);
                for jDime = 1:nDime
                    invJ_JI   = invJ(:,jDime,:,igaus);
                    dShape_KJ = dShapes(jDime,:);
                    dSDx_KI   = bsxfun(@times, invJ_JI,dShape_KJ);
                    dShapeDx(:,:,:,igaus) = dShapeDx(:,:,:,igaus) + dSDx_KI;
                end
            end
            dNdx = dShapeDx;
        end

        function gradFun = computeGradient(obj, quad)
            dNdx = obj.computeCartesianDerivatives(quad);
            nDimf = obj.ndimf;
            nDims = size(dNdx, 1); % derivX, derivY (mesh-related?)
            nNode = size(dNdx, 2);
            nElem = size(dNdx, 3);
            nGaus = size(dNdx, 4);
            
            grad = zeros(nDims,nDimf, nElem, nGaus);
            for iGaus = 1:nGaus
                dNdx_g = dNdx(:,:,:,iGaus);
                for iDims = 1:nDims
                    for iNode = 1:nNode
                        dNdx_i = squeeze(dNdx_g(iDims, iNode,:));
                        nodes = obj.connec(:,iNode);
                        f = obj.fValues(nodes,:);
                        p = (dNdx_i.*f)';
                        pp(1,:,:) = p;
                        grad(iDims,:,:,iGaus) = grad(iDims,:,:,iGaus) + pp;
                    end
                end
            end
            fVR = reshape(grad, [nDims*nDimf,nElem, nGaus]);
            s.fValues = permute(fVR, [1 3 2]);
%             s.ndimf      = nDimf;
            s.quadrature = quad;
            gradFun = FGaussDiscontinuousFunction(s);
        end

        function symGradFun = computeSymmetricGradient(obj,quad)
            grad = obj.computeGradient(quad);
            nDimf = obj.ndimf;
            nDims = size(grad.fValues, 1)/nDimf;
            nGaus = size(grad.fValues, 2);
            nElem = size(grad.fValues, 3);

            gradReshp = reshape(grad.fValues, [nDims,nDimf,nGaus,nElem]);
            gradT = permute(gradReshp, [2 1 3 4]);
            symGrad = 0.5*(gradReshp + gradT);
            
            s.fValues    = reshape(symGrad, [nDims*nDimf,nGaus,nElem]);
            s.quadrature = quad;
            symGradFun = FGaussDiscontinuousFunction(s);
        end

        function plot(obj) % 2D domains only
            if  strcmp(obj.order,'LINEAR')
                switch obj.mesh.type
                case {'TRIANGLE','QUAD'}
                    x = obj.coord(:,1);
                    y = obj.coord(:,2);
                    figure()
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
                    figure()
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
            if nargin == 2; software = 'GiD'; end
            sF.fValues = obj.fValues;
            sF.mesh = obj.mesh;
            p1 = P1Function(sF);
            s.mesh = obj.mesh;
            s.fun = {p1};
            s.type = software;
            s.filename = filename;
            p = FunctionPrinter.create(s);
            p.print();
        end
        
        function v = computeL2norm(obj)
            s.type     = 'ScalarProduct';
            s.quadType = 'QUADRATICMASS';
            s.mesh     = obj.mesh;
            int = Integrator.create(s);
            ff  = int.compute(obj,obj);
            v   = sqrt(ff);
        end

    end

    methods (Access = public, Static)

        function pL = create(mesh, ndimf, ord)
            s.mesh    = mesh;
            s.order   = ord;
            s.ndimf   = ndimf;
            c = DOFsComputer(s);
            c.computeDofs();
            s.fValues = zeros(c.getNumberDofs(),ndimf);
            pL = LagrangianFunction(s);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh = cParams.mesh;
            obj.fValues = cParams.fValues;
            obj.ndimf   = size(cParams.fValues,2);
            obj.order   = cParams.order;
        end

        function createInterpolation(obj)
            m.type = obj.mesh.type;
            obj.interpolation = Interpolation.create(m,obj.orderTextual());
        end

        function createDOFCoordConnec(obj)
            s.mesh          = obj.mesh;
            s.interpolation = obj.interpolation;
            s.order = obj.order;
            s.ndimf = obj.ndimf;
            c = DOFsComputer(s);
            c.computeDofs();
            c.computeCoord();
            obj.coord  = c.getCoord();
            obj.connec = c.getDofs();
        end
        
        function ord = orderTextual(obj)
            switch obj.order
                case 'P1'
                    ord = 'LINEAR';
                case 'P2'
                    ord = 'QUADRATIC';
                case 'P3'
                    ord = 'CUBIC';
            end
        end

    end

end