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
            % nNode  = obj.interpolation.nnode;
            % nGaus  = size(shapes,2);
            % nF     = size(obj.fValues,2);
            % nElem  = size(obj.dofConnec,1);
            % fxV = zeros(nF,nGaus,nElem);
            % nodes = obj.getDofConnecByVector();
            % for iGaus = 1:nGaus
            %     for iNode = 1:nNode
            %         node = nodes(:,iNode);
            %         Ni = shapes(iNode,iGaus);
            %         fi = obj.fValues(node,:);
            %         f(:,1,:) = Ni*fi';
            %         fxV(:,iGaus,:) = fxV(:,iGaus,:) + f;
            %     end
            % end     
           func = obj.getFvaluesDisc2();

          % fxV2 = fxV;

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

          %  norm(fxV2(:)-fxV(:))/norm(fxV2(:))

        end

        function fVals = getFvaluesDisc2(obj)
           nDimF  = size(obj.fValues,2);            
         %   fVals = reshape(obj.fValues',nDimF,[],obj.mesh.nelem);
            node = obj.getDofConnecByVector();         
            nnodeElem = size(node,2);
            for iDim = 1:nDimF
                fI = obj.fValues(:,iDim);
                for iNode = 1:nnodeElem
                    %iDof  = (iNode-1)*obj.ndimf+iDim;
%                    dofs  = obj.dofConnec(:,iDof);
                    dof   = node(:,iNode);

                    fVals(iDim,iNode,:) = fI(dof);
                            
                   % fVals(iDim,iNode,:) = fV;

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

        function div = computeDiv(obj)
            s.operation = @(xV) obj.computeDivFun(xV);
            s.ndimf     = 1;
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

        function plotVector(obj,varargin)
            if size(varargin,1) == 1
                n = varargin{1};
            else 
                n = 2;
            end
            figure();
            x = obj.mesh.coord(1:n:end,1);
            y = obj.mesh.coord(1:n:end,2);
            fX = obj.fValues(1:n:end,1);
            fY = obj.fValues(1:n:end,2);
            quiver(x, y, fX, fY, 'AutoScale', 'on', 'LineWidth', 1.5);              
            axis equal;  
            box on;     
            xlim([min(x), max(x)]);
            ylim([min(y), max(y)]);
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

        function f = copy(obj)
            s.ndimf = obj.ndimf;
            s.order = obj.order;
            s.mesh      = obj.mesh;
            s.fValues   = obj.fValues;
            s.dofConnec = obj.dofConnec;
            s.dofCoord  = obj.dofCoord;
            s.dofs.getNumberDofs = size(obj.dofCoord,1);
            f = LagrangianFunction(s);
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

        function r = mtimes(a,b)
            aOp = DomainFunction.computeOperation(a);
            bOp = DomainFunction.computeOperation(b);
            s.operation = @(xV) pagemtimes(aOp(xV),bOp(xV));
            r = DomainFunction(s);
        end

        function r = times(a,b)
            aOp = DomainFunction.computeOperation(a);
            bOp = DomainFunction.computeOperation(b);
            ndimfA = DomainFunction.computeFieldDimension(a);
            ndimfB = DomainFunction.computeFieldDimension(b);
            s.operation = @(xV) aOp(xV).*bOp(xV);
            s.ndimf = max(ndimfA,ndimfB);
            r = DomainFunction(s);
        end

        function f = power(f1,b)
            s.operation = @(xV) (f1.evaluate(xV)).^b;
            s.ndimf = f1.ndimf;
            f = DomainFunction(s);
        end

        function r = rdivide(a,b)
            aOp = DomainFunction.computeOperation(a);
            bOp = DomainFunction.computeOperation(b);
            s.operation = @(xV) aOp(xV)./bOp(xV);
            r = DomainFunction(s);
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
