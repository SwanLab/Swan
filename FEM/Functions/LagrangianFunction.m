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

        function node = getDofConnecByVector(obj)
          nNode = obj.interpolation.nnode;
          for iNode = 1:nNode
            iDof   = (iNode-1)*obj.ndimf+1;              
            node(:,iNode) = (obj.dofConnec(:,iDof)-1)/obj.ndimf+1;
          end
        end        

        function fxV = evaluate(obj, xV)
            shapes = obj.interpolation.computeShapeFunctions(xV);
            nNode  = obj.interpolation.nnode;
            nGaus  = size(shapes,2);
            nF     = size(obj.fValues,2);
            nElem  = size(obj.dofConnec,1);
            fxV = zeros(nF,nGaus,nElem);
            nodes = obj.getDofConnecByVector();
            for iGaus = 1:nGaus
                for iNode = 1:nNode
                    node = nodes(:,iNode);
                    Ni = shapes(iNode,iGaus);
                    fi = obj.fValues(node,:);
                    f(:,1,:) = Ni*fi';
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

        function dNdx  = evaluateCartesianDerivatives(obj,xV)
            if ~isequal(xV,obj.xVOld) || isempty(obj.dNdxOld)
                nElem = size(obj.dofConnec,1);
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
                obj.dNdxOld = dNdx;
                obj.xVOld   = xV;
            else
                dNdx = obj.dNdxOld;
            end
        end                

        function fVals = getFvaluesDisc(obj)
            nDimF  = size(obj.fValues,2);            
            fVals = reshape(obj.fValues',nDimF,[],obj.mesh.nelem);
            %obj.fValuesDisc = fVals;
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

        function dN = computeShapeDerivatives(obj, xV)
            dN = obj.interpolation.computeShapeDerivatives(xV);
        end

        function grad = computeGrad(obj)
            s.operation = @(xV) obj.computeGradFun(xV);
            s.ndimf     = obj.mesh.ndim*obj.ndimf;
            grad        = DomainFunction(s);
        end

        function curl = computeCurl(obj)
            s.operation = @(xV) obj.computeCurlFun(xV);
            s.ndimf     = 1;
            curl        = DomainFunction(s);            
        end

        function div = computeDiv(obj)
            s.operation = @(xV) obj.computeDivFun(xV);
            s.ndimf     = 1;
            div         = DomainFunction(s);                   
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
            ord = str2double(obj.order(end));
        end

        function plot(obj) % 2D domains only
            s.coord   = obj.getDofCoord();
            s.fValues = obj.fValues;
            s.ndimf   = obj.ndimf;
            switch obj.getOrderTextual(obj.order)
                case 'LINEAR'     
                    for iDim = 1:obj.ndimf
                        connecf{iDim} = obj.getDofConnecByVector();
                        coordf{iDim}  = obj.getDofCoordByVector(iDim);
                    end
                    s.connec = connecf;
                    s.coord  = coordf;
                    
                    lP = LagrangianPlotter(s);
                    lP.plot();                    
                case {'QUADRATIC','CUBIC'}
                    %better to remesh (now only plotting the linear part)                    
                    s.connec = obj.mesh.connec;                         
                    lP = LagrangianPlotter(s);
                    lP.plot();                    
                otherwise
                    f = obj.project('P1D');
                    f.plot()
            end

        end


        
        function cV = getDofCoordByVector(obj,dimf)
            cV = obj.getDofFieldByVector(dimf,obj.dofCoord);
        end

        function fV = getDofFieldByVector(obj,dimf,field)   
          ndimf = size(field,2);
          for iDim = 1:ndimf
                fieldD = field(:,iDim);
                fResh  = reshape(fieldD',obj.ndimf,[]);
                fV(:,iDim) = fResh(dimf,:);
          end         
        end        

        function fV = getfValuesAsVector(obj,dimf)    
            f = obj.fValues;
            fV = obj.getDofFieldByVector(dimf,f);
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

        function fFine = refine(obj,mFine) %Only for first order
            fNodes  = obj.fValues;
            fEdges  = obj.computeFunctionInEdges(obj.mesh, fNodes);
            fAll    = [fNodes;fEdges];
            s.mesh    = mFine;
            s.fValues = fAll;
            s.order   = obj.order;
            fFine = LagrangianFunction(s);
        end

        function f = copy(obj)
            f = obj.create(obj.mesh,obj.ndimf,obj.order);
            f.fValues = obj.fValues;
            f.setXvOld(obj.xVOld);
            f.setdNdxOld(obj.dNdxOld);            
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

        function s = times(f1,f2)
            s.operation = @(xV) f1.evaluate(xV) .* f2.evaluate(xV);
            s.ndimf = max(f1.ndimf,f2.ndimf);
            s = DomainFunction(s);
        end

        function f = power(f1,b)
            s.operation = @(xV) (f1.evaluate(xV)).^b;
            s.ndimf = f1.ndimf;
            f = DomainFunction(s);
        end

        function s = rdivide(f,b)
            res = copy(f);
            res.fValues = f.fValues ./ b;
            s = res;
        end

        function f = mrdivide(f1,f2)
            s.operation = @(xV) f1.evaluate(xV)./f2.evaluate(xV);
            s.ndimf = max(f1.ndimf,f2.ndimf);
            f = DomainFunction(s);            
        end

        function f = exp(f)
            s.operation = @(xV) exp(f.evaluate(xV));
            s.ndimf = f.ndimf;
            f = DomainFunction(s);
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
                case 'P1'
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
            fV    = obj.getValuesByElem();
            fV    = permute(fV,[1 2 4 3]);
            gradF = pagemtimes(dNdx,fV);
        end

        function divF = computeDivFun(obj,xV)
            nP = size(xV,2);
            dNdx  = obj.evaluateCartesianDerivatives(xV);
            fV    = obj.getValuesByElem();
            fV    = permute(fV,[1 2 4 3]);
            fV    = pagetranspose(fV);
            fV    = repmat(fV,[1 1 nP 1]);
            divF(1,:,:) = squeeze(bsxfun(@(A,B) sum(A.*B, [1 2]), fV,dNdx));        
        end

        function lapF = computeLaplacianFun(obj,xV)
            gradF = Grad(obj);
            gradF = gradF.project('P1',obj.mesh);
            lapF  = Divergence(gradF); 
            t.project('P1',obj.mesh).plot()              
        end

        function curlF = computeCurlFun(obj,xV)
            dNdx  = obj.evaluateCartesianDerivatives(xV);
            fV    = obj.getValuesByElem();
            fV    = obj.createOrthogonal(fV);
            fV    = permute(fV,[1 2 4 3]);
            gradF = pagemtimes(dNdx,fV);
            divF  = squeeze(sum(gradF,1));
            curlF  = divF;
        end       

        function fP = createOrthogonal(obj,f)
            fP = zeros(size(f));
            fP(:,1,:) = f(:,2,:);
            fP(:,2,:) = -f(:,1,:);
        end


       function fV = getValuesByElem(obj)
            connec = obj.getDofConnec();
            nNodeE = obj.interpolation.nnode;
            nElem  = obj.mesh.nelem;
            nDimf  = obj.ndimf;
            fV = zeros(nNodeE,nDimf,nElem);
            f  = reshape(obj.fValues', [1 obj.nDofs]);
            for iNode = 1:nNodeE
                for iDim = 1:nDimf
                    iDofE = nDimf*(iNode-1)+iDim;
                    dof = connec(:,iDofE);
                    fV(iNode,iDim,:) = f(dof);
                end
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
