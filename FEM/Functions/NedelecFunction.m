classdef NedelecFunction < FeFunction

    properties (GetAccess = public, SetAccess = private)
        nDofs
        nDofsElem
    end

    properties (Access = private)
        interpolation
        connec
    end

    methods (Access = public)

        function obj = NedelecFunction(cParams)
            obj.init(cParams);
            obj.createInterpolation();

            if not(contains(fieldnames(cParams),'dofs'))
                obj.createDOFCoordConnec();
            else
                obj.connec = cParams.dofs.getDofs();
                obj.nDofs = cParams.dofs.getNumberDofs();
            end
        end

        function fxV = evaluate(obj, xV)
            shapes = obj.interpolation.computeShapeFunctions(xV);
            nNode  = size(shapes,1);
            nGaus  = size(shapes,2);
            nDim   = obj.mesh.ndim;
            nF     = size(obj.fValues,2);
            nElem  = size(obj.connec,1);
            fxV = zeros(nDim,nGaus,nElem);

            shapesTestMapped  = obj.mapFunction(shapes, xV);

            for iGaus = 1:nGaus
                for iNode = 1:nNode
                    node = (obj.connec(:,(iNode-1)*obj.ndimf+1)-1)/obj.ndimf+1;
                    Ni   = reshape(squeeze(shapesTestMapped(iNode,iGaus,:,:))',nDim,[])';
                    fi   = obj.fValues(node,:);
                    f(1:nDim,1,:) = Ni'.*fi';
                    fxV(:,iGaus,:) = fxV(:,iGaus,:) + f;
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

        function mapF = mapFunction(obj, F, ~)
            mapF = zeros([size(F),obj.mesh.nelem]);
            J = obj.mesh.computeJacobian(0);

            JGlob = pagetranspose(pageinv(J));
            sides = obj.computeSidesOrientation();

            for idof= 1:obj.nDofsElem
                s(1,1,1,:) = sides(:,idof);
                mapF(idof,:,:,:,:) = pagemtimes(squeeze(F(idof,:,:,:)),JGlob).*s;
            end
        end

        function mapF = mapDerivFunction(obj, F, ~)
            mapF = zeros([size(F),obj.mesh.nelem]);
            J = obj.mesh.computeJacobian(0);

            JGlob = pagetranspose(pageinv(J));
            sides = obj.computeSidesOrientation();

            for idof= 1:obj.nDofsElem
                s(1,1,1,:) = sides(:,idof);
                mapF(idof,:,:,:,:) = pagemtimes(squeeze(F(idof,:,:,:)),JGlob).*s;
            end
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
                pl = NedelecPlotter();
                s.func = obj;
                s.mesh = obj.mesh;
                s.interpolation = obj.interpolation;
                pl.plot(s);
            end
        end

        % function dofConnec = computeDofConnectivity(obj)
        %     conne  = obj.connec;
        %     nDimf  = obj.ndimf;
        %     nNode  = size(conne, 2);
        %     nDofsE = nNode*nDimf;
        %     dofsElem  = zeros(nDofsE,size(conne,1));
        %     for iNode = 1:nNode
        %         for iUnkn = 1:nDimf
        %             idofElem   = nDimf*(iNode - 1) + iUnkn;
        %             globalNode = conne(:,iNode);
        %             idofGlobal = nDimf*(globalNode - 1) + iUnkn;
        %             dofsElem(idofElem,:) = idofGlobal;
        %         end
        %     end
        %     dofConnec = dofsElem;
        % end

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
            fFine = NedelecFunction(s);
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

        function ord = getOrderNum(obj)
            ord = 1; % NO
        end


        function loc = computeLocPointEdgeRef(obj)
            type = obj.mesh.type;
            switch type
                case 'LINE'
                    loc = [1 2];
                case 'TRIANGLE'
                    loc = [1 2 3];
                case 'QUAD'
                    loc = [1 2 3 4];
                case 'TETRAHEDRA'
                    loc = [1 1 1 2 2 3];
                case 'HEXAHEDRA'
                    loc = [1 4 1 2 2 3 3 4 5 8 6 7];
            end
        end


        function sides = computeSidesOrientation(obj)
            locPointEdge = squeeze(obj.mesh.edges.localNodeByEdgeByElem(:,:,1));
            sides = zeros(obj.mesh.nelem,obj.mesh.edges.nEdgeByElem);
            locPointEdgeRef = obj.computeLocPointEdgeRef();

            for iEdge = 1:obj.mesh.edges.nEdgeByElem
                sides(:,iEdge) = (locPointEdge(:,iEdge)==locPointEdgeRef(iEdge)).*2-1;
            end
        end
        
    end

    methods (Access = public, Static)

        function pL = create(mesh, ndimf, ord)
            s.mesh    = mesh;
            s.order   = ord;
            s.ndimf   = ndimf;
            c = DOFsComputerNedelec(s);
            c.computeDofs();
            s.fValues = zeros(c.getNumberDofs()/ndimf,ndimf);
            s.dofs = c;
            pL = NedelecFunction(s);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh    = cParams.mesh;
            obj.fValues = cParams.fValues;
            obj.ndimf   = size(cParams.fValues,2);
            % obj.order   = cParams.order;
        end

        function createInterpolation(obj)
            type = obj.mesh.type;
            obj.interpolation = Interpolation.create(type,'Nedelec');
            obj.nDofsElem = obj.ndimf*obj.interpolation.nnode;
        end

        function createDOFCoordConnec(obj)
            s.mesh          = obj.mesh;
            s.interpolation = obj.interpolation;
            s.ndimf         = obj.ndimf;
            c = DOFsComputerNedelec(s);
            c.computeDofs();
            obj.connec = c.getDofs();
            obj.nDofs = c.getNumberDofs();
        end

        function f = computeFunctionInEdges(~,m,fNodes)
            s.edgeMesh = m.computeEdgeMesh();
            s.fNodes   = fNodes;
            eF         = EdgeFunctionInterpolator(s);
            f = eF.compute();
        end

    end

end