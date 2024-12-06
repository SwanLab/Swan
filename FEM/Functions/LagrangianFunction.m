classdef LagrangianFunction < FeFunction

    properties (GetAccess = public, SetAccess = private)
        nDofs
        nDofsElem
    end

    properties (Access = private)
        interpolation
        dofCoord
        dofConnec

       dNdxOld
       xVOld     
    end

    methods (Access = public)

        function obj = LagrangianFunction(cParams)
            obj.init(cParams);
            obj.createInterpolation();

            if not(contains(fieldnames(cParams),'dofs'))
                obj.createDOFCoordConnec();
            else
                obj.dofConnec = cParams.dofConnec;
                obj.dofCoord  = cParams.dofCoord;
                obj.nDofs = cParams.dofs.getNumberDofs();
            end
        end

        function fxV = evaluate(obj, xV)
            shapes = obj.interpolation.computeShapeFunctions(xV);
            func   = obj.getFvaluesByElem();
            nNode  = size(shapes,1);
            nGaus  = size(shapes,2);
            nF     = size(func,1);
            nElem  = size(func,3);
            fxV = zeros(nF,nGaus,nElem);
            for kNode = 1:nNode
                shapeKJ = shapes(kNode,:,:);
                fKJ     = func(:,kNode,:);
                f = bsxfun(@times,shapeKJ,fKJ);
                fxV = fxV + f;
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

        function dNdx  = evaluateCartesianDerivatives(obj,xV)
            if ~isequal(xV,obj.xVOld) || isempty(obj.dNdxOld)
                nElem = size(obj.dofConnec,1);
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
                obj.dNdxOld = dNdx;
                obj.xVOld   = xV;
            else
                dNdx = obj.dNdxOld;
            end
        end         

       function fVals = getFvaluesByElem(obj)
            nDimF     = obj.ndimf;
            nNode     = obj.interpolation.nnode;
            nElem     = size(obj.mesh.connec, 1);            
            iDof      = (0:nNode-1)*obj.ndimf + 1;
            node      = (obj.dofConnec(:, iDof) - 1) / nDimF + 1;
            fAll      = obj.fValues(node(:), :);
            fReshaped = reshape(fAll, nElem, nNode, nDimF);
            fVals     = permute(fReshaped, [3, 2, 1]);            
        end           

        function c = getDofCoord(obj)
            c = obj.dofCoord;
        end

        function c = getDofConnec(obj)
            c = obj.dofConnec;
        end

        function N = computeShapeFunctions(obj, xV)
            N = obj.interpolation.computeShapeFunctions(xV);
        end

        function grad = computeGrad(obj)
            s.operation = @(xV) obj.computeGradFun(xV);
            s.ndimf     = obj.mesh.ndim*obj.ndimf;
            s.mesh      = obj.mesh;
            grad        = DomainFunction(s);
        end

        function div = computeDiv(obj)
            s.operation = @(xV) obj.computeDivFun(xV);
            s.ndimf     = 1;
            s.mesh      = obj.mesh;
            div         = DomainFunction(s);                   
        end        

        function curl = computeCurl(obj) %only for 2D
            fOrth = obj.createOrthogonalVector();            
            curl  = Divergence(fOrth);            
        end

        function setdNdxOld(obj,dNdx)
            obj.dNdxOld = dNdx;
        end

        function setXvOld(obj,xV)
            obj.xVOld = xV;
        end           
        
        function ord = orderTextual(obj)
            ord = obj.getOrderTextual(obj.order);
        end

        function ord = getOrderNum(obj)
            switch obj.order
                case 'P1D'
                    ord = 1;
                otherwise
                ord = str2double(obj.order(end));
            end
        end

        function plot(obj) % 2D domains only
            switch obj.getOrderTextual(obj.order)
                case 'LINEAR'
                    figure()
                    connecP = obj.getDofConnecByVector();
                    for iDim = 1:obj.ndimf                                                
                        subplot(1,obj.ndimf,iDim);
                        coordP  = obj.getDofFieldByVector(iDim,obj.dofCoord);                                                
                        x  = coordP(:,1);
                        y  = coordP(:,2);
                        z  = double(obj.fValues(:,iDim));
                        a = trisurf(connecP,x,y,z);
                        view(0,90)
                        %colorbar
                        shading interp
                        a.EdgeColor = [0 0 0];
                        title(['dim = ', num2str(iDim)]);
                    end
                otherwise
                    f = obj.project('P1D');
                    f.plot()
            end
        end

        function plotVector(obj,varargin) %only for linear
            if size(varargin, 1) == 1, n = varargin{1}; else, n = 2; end
            figure();
            coord = obj.getDofFieldByVector(1,obj.dofCoord); 
            x = coord(1:n:end,1);
            y = coord(1:n:end,2);
            fX = obj.fValues(1:n:end,1);
            fY = obj.fValues(1:n:end,2);
            quiver(x, y, fX, fY, 'AutoScale', 'on', 'LineWidth', 1.5);              
            axis equal;  
            box on;     
            xlim([min(x), max(x)]);
            ylim([min(y), max(y)]);
        end

        function fV = getDofFieldByVector(obj,dimf,field)   
          ndimf = size(field,2);
          for iDim = 1:ndimf
                fieldD = field(:,iDim);
                fResh  = reshape(fieldD',obj.ndimf,[]);
                fV(:,iDim) = fResh(dimf,:);
          end         
        end        

        function dof = getDofsFromCondition(obj, condition)
            nodes = condition(obj.dofCoord);
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
            int = Integrator.compute(obj.*obj,obj.mesh,2);
            v   = sqrt(int);
        end

        function f = createOrthogonalVector(obj) %only in 2D and vector
            f = obj.copy();
            f.nDofs = obj.nDofs;
            f.fValues(:,1) = obj.fValues(:,2);
            f.fValues(:,2) = -obj.fValues(:,1);
        end

        function fFine = refine(obj,mFine)%Only for first order
            switch obj.order
                case 'P1'
                    fNodes  = obj.fValues;
                    fEdges  = obj.computeFunctionInEdges(obj.mesh, fNodes);
                    fAll    = [fNodes;fEdges];
                    s.mesh    = mFine;
                    s.fValues = fAll;
                    s.order   = obj.order;
                    fFine = LagrangianFunction(s);
                case 'P1D'
                    P1Dref = P1Refiner(obj,mFine);
                    fFine  = P1Dref.compute();
            end
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

        function s = plus(a,b)
            if isa(a, 'LagrangianFunction')
                res = copy(a);
                val1 = a.fValues;
            else
                val1 = a;
            end
            if isa(b, 'LagrangianFunction')
                res = copy(b);
                val2 = b.fValues;
            else
                val2 = b;
            end

            res.fValues = val1 + val2;
            s = res;
        end

        function s = minus(a,b)
            if isa(a, 'LagrangianFunction')
                res = copy(a);
                val1 = a.fValues;
            else
                val1 = a;
            end
            if isa(b, 'LagrangianFunction')
                res = copy(b);
                val2 = b.fValues;
            else
                val2 = b;
            end
            res.fValues = val1 - val2;
            s = res;
        end

        

    end

    methods (Access = public, Static)

        function pL = create(mesh, ndimf, ord)
            s.mesh    = mesh;
            s.order   = ord;
            s.ndimf   = ndimf;
            s.interpolation = Interpolation.create(mesh.type,LagrangianFunction.getOrderTextual(ord));
            dofs = DOFsComputer(s);
            dofs.computeDofs();
            dofs.computeCoord();            
            s.dofCoord  = dofs.getCoord();
            s.dofConnec = dofs.getDofs();
            s.fValues = zeros(dofs.getNumberDofs()/ndimf,ndimf);            
            pL = LagrangianFunction(s);
        end
        
        function ord = getOrderTextual(order)
            switch order
                case 'P0'
                    ord = 'CONSTANT';
                case {'P1','P1D'}
                    ord = 'LINEAR';
                case 'P2'
                    ord = 'QUADRATIC';
                case 'P3'
                    ord = 'CUBIC';
            end        
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
            obj.interpolation = Interpolation.create(type,obj.getOrderTextual(obj.order));
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
            obj.dofCoord  = c.getCoord();
            obj.dofConnec = c.getDofs();
            obj.nDofs  = c.getNumberDofs();
        end

        function gradF = computeGradFun(obj,xV)
            dNdx  = obj.evaluateCartesianDerivatives(xV);
            fV    = obj.getFvaluesByElem();            
            fV    = permute(fV,[2 1 4 3]);
            gradF = pagemtimes(dNdx,fV);
        end

        function divF = computeDivFun(obj,xV)
            nP = size(xV,2);
            dNdx  = obj.evaluateCartesianDerivatives(xV);
            fV    = obj.getFvaluesByElem(); 
            fV    = permute(fV,[2 1 4 3]);
            fV    = pagetranspose(fV);
            fV    = repmat(fV,[1 1 nP 1]);
            divF(1,:,:) = squeeze(bsxfun(@(A,B) sum(A.*B, [1 2]), fV,dNdx));        
        end

        function lapF = computeLaplacianFun(obj,xV)
            gradF = Grad(obj);
            gradF = gradF.project('P1',obj.mesh);
            lapF  = Divergence(gradF); 
        end

        function node = getDofConnecByVector(obj)
            nNode = obj.interpolation.nnode;
            nElem = size(obj.mesh.connec,1);
            node  = zeros(nElem,nNode);
            for iNode = 1:nNode
                iDof   = (iNode-1)*obj.ndimf+1;
                node(:,iNode) = (obj.dofConnec(:,iDof)-1)/obj.ndimf+1;
            end
        end

      

        function f = computeFunctionInEdges(obj,m,fNodes)
            s.edgeMesh = m.computeEdgeMesh();
            s.fNodes   = fNodes;
            eF         = EdgeFunctionInterpolator(s);
            f = eF.compute()';
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
