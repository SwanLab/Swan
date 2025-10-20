classdef LagrangianFunction < FeFunction

    properties (GetAccess = public, SetAccess = private)
        nDofs
        nDofsElem
        dofCoord
        dofConnec
        % bFun
    end

    properties (Access = private)
        interpolation
        dNdxOld
        xVOlddN
    end

    methods (Access = public, Static)

        function fun = create(mesh,ndimf,ord)
            s.mesh    = mesh;
            s.order   = ord;
            s.ndimf   = ndimf;
            fun = LagrangianFunction(s);
        end

        function ord = getOrderTextual(order)
            switch order
                case 'P0';         ord = 'CONSTANT';
                case {'P1','P1D'}; ord = 'LINEAR';
                case 'P2';         ord = 'QUADRATIC';
                case 'P3';         ord = 'CUBIC';
            end        
        end

    end

    methods (Access = public)

        function obj = LagrangianFunction(cParams)
            obj.init(cParams);
            obj.createInterpolation();
            obj.createDOFCoordConnec();
            obj.createFValues(cParams);
        end

        function setFValues(obj,fV)
            obj.fValues = fV;
            obj.fxVOld  = [];
        end

        function fVElem = getFvaluesByElem(obj)
            nDimTotal = obj.ndimfTotal;
            nNode     = obj.interpolation.nnode;
            nElem     = obj.mesh.nelem;
            fVdofs    = obj.fValues(obj.dofConnec)';
            fVElem = reshape(fVdofs,nDimTotal,nNode,nElem);
        end

        function fVNode = getFValuesByNode(obj)
            nnodes  = obj.mesh.nnodes;
            fVNode = reshape(obj.fValues,obj.ndimfTotal,nnodes)';
        end
        
       %  function fxV = sampleFunction(obj,xP,cells)
       %      shapes  = obj.interpolation.computeShapeFunctions(xP);
       %      nNode   = size(shapes,1);
       %      nF      = size(obj.fValues,2);
       %      nPoints = size(xP,2);
       %      fxV = zeros(nF,nPoints);
       %      for iF = 1:nF
       %          for iNode = 1:nNode
       %              node = obj.mesh.connec(cells,iNode);
       %              Ni = shapes(iNode,:)';
       %              fi = obj.fValues(node,:);
       %              f(1,:) = fi.*Ni;
       %              fxV(iF,:) = fxV(iF,:) + f;
       %          end
       %      end
       %  end
       % 

       function N = computeShapeFunctions(obj, xV)
            N = obj.interpolation.computeShapeFunctions(xV);
       end

        function dNdx  = evaluateCartesianDerivatives(obj,xV)
            if ~isequal(xV,obj.xVOlddN) || isempty(obj.dNdxOld)
                invJ    = Inv(Jacobian(obj.mesh));
                N       = obj.interpolation;
                dShapes = invJ*obj.Gradient(N);
                dNdx    = dShapes.evaluate(xV);
                obj.dNdxOld = dNdx;
                obj.xVOlddN = xV;
            else
                dNdx = obj.dNdxOld;
            end
        end 

       function gradF = computeGrad(obj,xV)
            dNdx  = obj.evaluateCartesianDerivatives(xV);
            fV = obj.getFvaluesByElem(); 
            fV = permute(fV,[2 1 4 3]);
            gradF = squeezeParticular(pagemtimes(dNdx,fV),[1 2]);
        end

        function div = computeDiv(obj)
            s.operation = @(xV) obj.computeDivFun(xV);
            s.ndimf     = 1;
            s.mesh      = obj.mesh;
            div         = DomainFunction(s);                   
        end        

        % function curl = computeCurl(obj) %only for 2D
        %     fOrth = obj.createOrthogonalVector();            
        %     curl  = Divergence(fOrth);            
        % end
       % 
       %  function setdNdxOld(obj,dNdx)
       %      obj.dNdxOld = dNdx;
       %  end
       % 
       %  function setXvdNOld(obj,xV)
       %      obj.xVOlddN = xV;
       %  end  
       % 
       %  function setXvfVOld(obj,xV)
       %      obj.xVOldfV = xV;
       %  end   
       % 
       %  function ord = orderTextual(obj)
       %      ord = obj.getOrderTextual(obj.order);
       %  end
       % 
       %  function ord = getOrderNum(obj)
       %      switch obj.order
       %          case 'P1D'
       %              ord = 1;
       %          otherwise
       %          ord = str2double(obj.order(end));
       %      end
       %  end
       % 
       %  function plot(obj) % 2D domains only
       %      plotFun = @(tri,x,y,z,iDim) obj.plotF(tri,x,y,z,iDim);
       %      obj.generalPlot(plotFun)
       %  end
       % 
       %  function plotContour(obj,varargin) % 2D domains only
       %  if size(varargin, 1) == 1, nC = varargin{1}; else, nC = 30; end            
       %     plotFun  = @(tri,x,y,z,iDim) obj.plotContourF(tri,x,y,z,iDim,nC);
       %     plotC    = @() obj.plotContour();
       %     obj.generalPlot(plotFun,plotC)
       %  end     
       % 
       %  function plotIsoLines(obj,varargin) % 2D domains only 
       %      if size(varargin, 1) == 1, nC = varargin{1}; else, nC = 16; end  
       %      figure()
       %      connecP  = obj.getDofConnecByVector();
       %      for iDim = 1:obj.ndimf
       %          coordP  = obj.getDofFieldByVector(iDim,obj.dofCoord);
       %          x  = coordP(:,1);
       %          y  = coordP(:,2);
       %          z  = double(obj.fValues(:,iDim));
       %          [~,h] = tricontour(connecP,x,y,z,nC);
       %          view(0,90)            
       %          set(h,'LineWidth',2);
       %          set(h, 'EdgeColor', 'k');
       %          view(0,90)
       %          hold on
       %      end
       %  end
       % 
       %  function plotVector(obj,varargin) %only for linear
       %      if size(varargin, 1) == 1, n = varargin{1}; else, n = 2; end
       %      figure();
       %      coord = obj.getDofFieldByVector(1,obj.dofCoord); 
       %      x = coord(1:n:end,1);
       %      y = coord(1:n:end,2);
       %      fX = obj.fValues(1:n:end,1);
       %      fY = obj.fValues(1:n:end,2);
       %      quiver(x, y, fX, fY, 'AutoScale', 'on', 'LineWidth', 1.5);              
       %      axis equal;  
       %      box on;     
       %      xlim([min(x), max(x)]);
       %      ylim([min(y), max(y)]);            
       %  end
       % 
       %  function fV = getDofFieldByVector(obj,dimf,field)   
       %    ndimf = size(field,2);
       %    for iDim = 1:ndimf
       %          fieldD = field(:,iDim);
       %          fResh  = reshape(fieldD',obj.ndimf,[]);
       %          fV(:,iDim) = fResh(dimf,:);
       %    end         
       %  end        
       % 
       %  function dof = getDofsFromCondition(obj, condition)
       %      nodes = condition(obj.dofCoord);
       %      iNode = find(nodes==1);
       %      dofElem = repmat(1:obj.ndimf, [length(iNode) 1]);
       %      dofMat = obj.ndimf*(iNode - 1) + dofElem;
       %      dof = sort(dofMat(:));
       %  end
       % 
       %  function print(obj, filename, software)
       %      if nargin == 2; software = 'Paraview'; end
       %      %             sF.fValues = obj.fValues;
       %      %             sF.mesh = obj.mesh;
       %      %             p1 = P1Function(sF);
       %      s.mesh = obj.mesh;
       %      s.fun = {obj};
       %      s.type = software;
       %      s.filename = filename;
       %      p = FunctionPrinter.create(s);
       %      p.print();
       %  end
       % 
       %  function [res, pformat] = getDataToPrint(obj)
       %      switch obj.order
       %          case 'P0'
       %              q = Quadrature.set(obj.mesh.type);
       %              q.computeQuadrature('LINEAR');
       %              nElem = size(obj.mesh.connec, 1);
       %              nGaus = q.ngaus;
       % 
       %              s.nDimf   = obj.ndimf;
       %              s.nData   = nElem*nGaus;
       %              s.nGroup  = nElem;
       %              s.fValues = obj.getFormattedP0FValues();
       %              fps = FunctionPrintingSettings(s);
       %              [res, pformat] = fps.getDataToPrint();
       % 
       %          otherwise
       %              nNods = size(obj.fValues, 1);
       %              s.nDimf   = obj.ndimf;
       %              s.nData   = nNods;
       %              s.nGroup  = nNods;
       %              s.fValues = obj.fValues;
       %              fps = FunctionPrintingSettings(s);
       %              [res, pformat] = fps.getDataToPrint();
       %      end
       %  end
       % 
       %  function f = createOrthogonalVector(obj) %only in 2D and vector
       %      f = copy(obj);
       %      f.fxVOld = [];
       %      f.nDofs = obj.nDofs;
       %      f.fValues(:,1) = obj.fValues(:,2);
       %      f.fValues(:,2) = -obj.fValues(:,1);
       %  end
       % 
       %  function f = getVectorFields(obj)
       %      for iDim = 1:obj.ndimf
       %          s.ndimf     = obj.ndimf;
       %          s.dofConnec = obj.getDofConnecByVector;
       %          s.dofCoord  = obj.getDofFieldByVector(iDim,obj.dofCoord);                                                
       %          s.fValues   = obj.fValues(:,iDim);
       %          s.mesh      = obj.mesh;
       %          s.order     = obj.order;
       %          f{iDim} = LagrangianFunction(s);                
       %      end
       % 
       %  end
       % 
       %  function fFine = refine(obj,mFine) %Only for first order
       %      switch obj.order
       %          case 'P1'
       %              fNodes  = obj.fValues;
       %              fEdges  = obj.computeFunctionInEdges(obj.mesh, fNodes);
       %              fAll    = [fNodes;fEdges];
       %              s.mesh    = mFine;
       %              s.fValues = fAll;
       %              s.order   = obj.order;
       %              fFine = LagrangianFunction(s);
       %          case 'P1D'
       %              P1Dref = P1Refiner(obj,mFine);
       %              fFine  = P1Dref.compute();
       %      end
       %  end
       % 
       %  function f = normalize(obj,type,varargin)
       %      fNorm = Norm(obj,type,varargin{:});
       %      f = obj.create(obj.mesh,obj.ndimf,obj.order);
       %      f.fValues = obj.fValues/fNorm;
       %  end
       % 
       %  function bF = restrictBaseToBoundary(obj,bMesh)
       %      if isempty(obj.bFun)
       %          obj.bFun = LagrangianFunction.create(bMesh,obj.ndimf,obj.order);
       %      end
       %      bF = obj.bFun;
       %  end
       % 
       %  % Operator overload
       % 
       %  function s = plus(a,b)
       %      if isnumeric(a) || isnumeric(b)
       %          if isa(a, 'LagrangianFunction')
       %              res = copy(a);
       %              val1 = a.fValues;
       %              fEv1 = a.fxVOld;
       %          else
       %              val1 = a;
       %              fEv1 = a;
       %          end
       %          if isa(b, 'LagrangianFunction')
       %              res = copy(b);
       %              val2 = b.fValues;
       %              fEv2 = b.fxVOld;
       %          else
       %              val2 = b;
       %              fEv2 = b;
       %          end
       %          if ~isempty(fEv1) && ~isempty(fEv2)
       %              res.fxVOld = fEv1 + fEv2;
       %          else
       %              res.fxVOld = [];
       %          end
       %          res.fValues = val1 + val2;
       %          s = res;
       %      else % a will be lagrangian, otherwise won't enter here              
       %          if isa(b, 'LagrangianFunction') & b.order == a.order
       %              res = copy(a);
       %              val1 = a.fValues;
       %              fEv1 = a.fxVOld;
       %              val2 = b.fValues;
       %              fEv2 = b.fxVOld;
       %              if ~isempty(fEv1) && ~isempty(fEv2)
       %                  res.fxVOld = fEv1 + fEv2;
       %              else
       %                  res.fxVOld = [];
       %              end
       %              res.fValues = val1 + val2;
       %              s = res;
       %          elseif isa(b, 'BaseFunction')
       %              s = plus@BaseFunction(a,b);
       %          end
       %      end
       %  end
       % 
       %  function s = minus(a,b)
       %      s = plus(a,-b);
       %  end
       % 
       %  function s = uminus(a)
       %      s = copy(a);
       %      s.setFValues(-a.fValues);
       %  end

    end

    methods (Access = protected)

        function fxV = evaluateNew(obj, xV)
            nGauss = size(xV,2);
            nElem  = obj.mesh.nelem;
            shapes = obj.interpolation.computeShapeFunctions(xV);
            fVal   = pagetranspose(obj.getFvaluesByElem());
            fxVvec = pagetranspose(tensorprod(shapes,fVal,1,1));
            fxVmat = reshape(fxVvec,[obj.ndimf,nGauss,nElem]);
            transposedDims = length(obj.ndimf):-1:1;
            extraDims      = length(obj.ndimf)+(1:2);
            fxV = permute(fxVmat,[transposedDims extraDims]);
        end        

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh    = cParams.mesh;
            obj.ndimf   = cParams.ndimf;
            obj.order   = cParams.order;
            obj.ndimfTotal = prod(obj.ndimf);
        end

        function createInterpolation(obj)
            type = obj.mesh.type;
            obj.interpolation = Interpolation.create(type,obj.getOrderTextual(obj.order));
            obj.nDofsElem = obj.ndimfTotal*obj.interpolation.nnode;
        end
        
        function createDOFCoordConnec(obj)
            s.mesh          = obj.mesh;
            s.interpolation = obj.interpolation;
            s.order         = obj.order;
            s.ndimf         = obj.ndimfTotal;
            c = DOFsComputer(s);
            c.computeDofs();
            c.computeCoord();
            obj.dofCoord  = c.getCoord();
            obj.dofConnec = c.getDofs();
            obj.nDofs  = c.getNumberDofs();
        end

        function createFValues(obj,cParams)
            obj.fValues = zeros(obj.nDofs,1);
            if isfield(cParams,'fValues')
                obj.fValues = cParams.fValues;
            end
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

        function lapF = computeLaplacianFun(obj,xV) %% (This should be a Function -> Laplacian) + (Not used anywhere)
            gradF = Grad(obj);
            gradF = gradF.project('P1',obj.mesh);
            lapF  = Divergence(gradF); 
        end
         
        % function node = getDofConnecByVector(obj)
        %     nNode = obj.interpolation.nnode;
        %     nElem = size(obj.mesh.connec,1);
        %     node  = zeros(nElem,nNode);
        %     for iNode = 1:nNode
        %         iDof   = (iNode-1)*obj.ndimf+1;
        %         node(:,iNode) = (obj.dofConnec(:,iDof)-1)/obj.ndimf+1;
        %     end
        % end
        % 
        % 
        % 
        % function f = computeFunctionInEdges(obj,m,fNodes)
        %     s.edgeMesh = m.computeEdgeMesh();
        %     s.fNodes   = fNodes;
        %     eF         = EdgeFunctionInterpolator(s);
        %     f = eF.compute()';
        % end
        % 
        % function fM = getFormattedP0FValues(obj)
        %     q = Quadrature.set(obj.mesh.type);
        %     q.computeQuadrature('LINEAR');
        %     fV = obj.evaluate(q.posgp);
        %     nGaus   = q.ngaus;
        %     nComp   = obj.ndimf;
        %     nElem   = size(obj.mesh.connec, 1);
        %     fM  = zeros(nGaus*nElem,nComp);
        %     for iStre = 1:nComp
        %         for iGaus = 1:nGaus
        %             rows = linspace(iGaus,(nElem - 1)*nGaus + iGaus,nElem);
        %             fM(rows,iStre) = fV(iStre,iGaus,:);
        %         end
        %     end
        % end
        % 
        % function plotF(obj,tri,x,y,z,iDim)
        %     a = trisurf(tri,x,y,z);
        %     view(0,90)
        %     %colorbar
        %     shading interp
        %     a.EdgeColor = [0 0 0];
        %     title(['dim = ', num2str(iDim)]);
        % end
        % 
        % function plotContourF(obj,tri,x,y,z,iDim,nC) % 2D domains only                       
        %     [c,h] = tricontour(tri,x,y,z,nC);
        %     view(0,90)            
        %     set(h,'LineWidth',5);
        %     view(0,90)
        %     %colorbar
        %     title(['dim = ', num2str(iDim)]);            
        %     c = reshape(c,2,4,[]);   
        %     nP = size(c,3);
        %     c2 = c(:,:,(2*nP/nC):end);
        %     c3 = reshape(c2,2,[]);
        %     %clabel([],h,'LabelSpacing',72,'Color','b','FontWeight','bold');            
        %     clabel(c3);            
        % end      
        % 
        % 
        % 
        % function generalPlot(obj,plotFun)
        %     switch obj.getOrderTextual(obj.order)
        %         case 'LINEAR'
        %             figure()
        %             connecP  = obj.getDofConnecByVector();
        %             for iDim = 1:obj.ndimf                                                
        %                 subplot(1,obj.ndimf,iDim);
        %                 coordP  = obj.getDofFieldByVector(iDim,obj.dofCoord);                                                
        %                 x  = coordP(:,1);
        %                 y  = coordP(:,2);
        %                 z  = double(obj.fValues(:,iDim));
        %                 plotFun(connecP,x,y,z,iDim);
        %             end
        %         otherwise
        %             f = obj.project('P1D');
        %             plot(f);
        %     end
        % end
        %
        function gradN = Gradient(obj,N) % This should be included in interpolation
            s.operation = @(xV) N.computeShapeDerivatives(xV);
            s.mesh      = obj.mesh;
            gradN       = DomainFunction(s);
        end

    end

end
