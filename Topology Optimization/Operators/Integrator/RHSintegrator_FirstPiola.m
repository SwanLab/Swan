classdef RHSintegrator_FirstPiola < RHSintegrator
    
    % \int Piola : Fdot

    properties (Access = private)
        mu, lambda
    end

    methods (Access = public)
        
        function obj = RHSintegrator_FirstPiola(cParams)
            obj.init(cParams)
            obj.setQuadratureOrder(cParams);
            obj.createQuadrature();
        end
        
        function rhs = compute(obj, fun, test)
            rhsElem = obj.computeElementalRHS(fun,test);
            rhs = obj.assembleIntegrand(rhsElem,test);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj, cParams)
            obj.mesh = cParams.mesh;
            obj.lambda = 1;
            obj.mu = 1;
            obj.quadratureOrder = 'LINEAR';
        end
        
        function rhsC = computeElementalRHS(obj, uFun, test)
            % Matrices
            quad = obj.quadrature;
            xG = quad.posgp;
            nPoints  = obj.quadrature.ngaus;
            nElem = obj.mesh.nelem;
            nDimG = obj.mesh.ndim;
            nDimf = uFun.ndimf;
            GradU = reshape(Grad(uFun).evaluate(xG),[nDimG,nDimf,nPoints, nElem]);
            I33 = zeros(size(GradU));
            I33(1,1,:,:) = 1;
            I33(2,2,:,:) = 1;

            I33(3,3,:,:) = 1;
            F = I33 + GradU; % deformation gradient
            invF = MatrixVectorizedInverter.computeInverse(F);
            invFt = permute(invF, [2 1 3 4]);
            jac(1,1,:,:)  = MatrixVectorizedInverter.computeDeterminant(F);
            
            piola = obj.mu*(F-invFt) + obj.lambda*(jac-1).*jac.*invFt;
            % can i subtract 1 to the jacobian directly?

            dV(1,1,:,:) = obj.mesh.computeDvolume(obj.quadrature);
            dNdx = test.evaluateCartesianDerivatives(obj.quadrature.posgp);
            nNodeE = size(dNdx,2);
            % piola:dNdx
            % piola: nDimG*nDimF*nGaus*nElem
            % dNdx:  nDimG*nNodE*nGaus*nElem
            % do we need to transpose piola to make it consistent?
            % ndimF*nDimG x nDimG*nNodeE

            piola = permute(piola, [2 1 3 4]);
            mult = pagemtimes(piola, dNdx);
            intI = mult.*dV;
            rhsC = squeezeParticular(sum(intI,3),3);
            rhsC = reshape(rhsC, [nDimf*nNodeE, nElem]);
            % rhsC = zeros(nNode*nDim,nElem);
        end

        function f = assembleIntegrand(obj, rhsElem, test)
            integrand = pagetranspose(rhsElem);
            connec = test.getConnec();
            nDofs = max(max(connec));
            nDofElem  = size(connec,2);
            f = zeros(nDofs,1);
            for idof = 1:nDofElem
                int = integrand(:,idof);
                con = connec(:,idof);
                f = f + accumarray(con,int,[nDofs,1],@sum,0);
            end
        end
        
    end
    
end
