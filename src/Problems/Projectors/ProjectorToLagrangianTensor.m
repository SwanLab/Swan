classdef ProjectorToLagrangianTensor < Projector
    
    properties (Access = private)
        order
    end
    
    methods (Access = public)

        function obj = ProjectorToLagrangianTensor(cParams)
            obj.order = cParams.projectorType;
        end

        function xFun = project(obj, x)
            LHS = obj.computeLHS(x);
            RHS = obj.computeRHS(x);
            xProj = LHS\RHS;
            s.mesh    = x.mesh;
            s.fValues = xProj;
            s.order = obj.order;
            xFun = LagrangianFunction(s);
        end

    end

    methods (Access = private)
        
        function LHS = computeLHS(obj,fun)
            switch obj.order
                case 'P0'
                    quad = Quadrature.create(fun.mesh,1);
                    dv = fun.mesh.computeDvolume(quad);
                    a = sum(dv(1,:),1);
                    a = repmat(a,1,1);
                    LHS = spdiags(a',0,length(a),length(a));
                otherwise
                    s.mesh  = fun.mesh;
                    s.test  = LagrangianFunction.create(fun.mesh, 1, obj.order);
                    s.trial = LagrangianFunction.create(fun.mesh, 1, obj.order);
                    s.type  = 'MassMatrix';
                    lhs = LHSIntegrator.create(s);
                    LHS = lhs.compute();
            end
        end

        function RHS = computeRHS(obj,fun)
            ord   = obj.createRHSQuadrature(fun);
            quad  = Quadrature.create(fun.mesh,ord);
            xV    = quad.posgp;
            fG    = fun.evaluate(xV);
            dV    = fun.mesh.computeDvolume(quad);
            test   = LagrangianFunction.create(fun.mesh,1,obj.order);
            N     = obj.computeShapeFunctionsBase(fun.mesh,test,xV);
            nDof  = size(N,3);
            nGaus = size(N,4);
            nElem = fun.mesh.nelem;
            int   = zeros(nDof,nElem);
            for iDof = 1:nDof
                for iGaus = 1:nGaus
                    dVg(:,1) = dV(iGaus, :);
                    fV   = squeezeParticular(fG(:,:,iGaus,:),3);
                    Ni   = N(:,:,iDof,iGaus);
                    Ni   = repmat(Ni,[1 1 nElem]);
                    fVNi = pagetensorprod(fV,Ni,[1 2],[1 2],2,2);
                    fNdV(1,:) = fVNi.*dVg;
                    int(iDof,:) = int(iDof,:) + fNdV;
                end
            end
            RHS = obj.assembleIntegrand(test,int);
        end

        function f = assembleIntegrand(obj,test,int)
            nDim     = test.mesh.ndim;
            nDofsEl  = test.nDofsElem;
            nElem    = test.mesh.nelem;
            int      = reshape(int,[nDim^2 nDofsEl nElem]);
            connec   = test.getDofConnec();
            ndofs    = max(max(connec));
            nDofElem = size(connec,2);
            f = zeros(ndofs,size(int,1));
            for iDof = 1:size(int,1)
                for iDofEl = 1:nDofElem
                    inti = squeeze(int(iDof,iDofEl,:));
                    con = connec(:,iDofEl);
                    f(:,iDof) = f(:,iDof) + accumarray(con,inti,[ndofs,1],@sum,0);
                end
            end
        end

        function N = computeShapeFunctionsBase(obj,mesh,test,xV)
            nDim    = mesh.ndim;
            nDofsEl = test.nDofsElem;
            nDof    = nDim^2*nDofsEl;
            nGauss  = size(xV,2);
            NxV     = test.computeShapeFunctions(xV);
            N       = zeros(nDim,nDim,nDof,nGauss);
            for i = 1:nDofsEl
                Ng = NxV(i,:);
                for j = 1:nDim
                    for k = 1:nDim
                        iDof = k + (j-1)*nDim + (i-1)*nDim^2;
                        N(j,k,iDof,:) = Ng;
                    end
                end
            end
        end

        function ord = createRHSQuadrature(obj, fun)
            ord = obj.determineQuadratureOrder(fun);
        end

    end

end