classdef LHSintegrator_SecondPiola < RHSintegrator
    
    % \int Piola : Fdot

    properties (Access = private)
        mu, lambda
    end

    methods (Access = public)
        
        function obj = LHSintegrator_SecondPiola(cParams)
            obj.init(cParams)
            obj.setQuadratureOrder(cParams);
            obj.createQuadrature();
        end
        
        function rhs = compute(obj, fun, test)
            rhsElem = obj.computeElementalLHS(fun,test);
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
        
        function lhsC = computeElementalLHS(obj, uFun, test)
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
            Ft = permute(F, [2 1 3 4]);
            C = pagemtimes(Ft,F);
            invC = MatrixVectorizedInverter.computeInverse(C);
            jac(1,1,:,:)  = MatrixVectorizedInverter.computeDeterminant(F);
            
            secPiola = obj.mu*(I33-invC) + obj.lambda*(jac-1).*jac.*invC;
            % can i subtract 1 to the jacobian directly?

            dV(1,1,:,:) = obj.mesh.computeDvolume(obj.quadrature);
            dNdx = test.evaluateCartesianDerivatives(obj.quadrature.posgp);
            nNodeE = size(dNdx,2);
            % piola:dNdx
            % piola: nDimG*nDimF*nGaus*nElem
            % dNdx:  nDimG*nNodE*nGaus*nElem
            % do we need to transpose piola to make it consistent?
            % ndimF*nDimG x nDimG*nNodeE

            dNdxT = permute(dNdx, [2 1 3 4]);
            secPiola = permute(secPiola, [2 1 3 4]);
            mult = pagemtimes(secPiola, dNdx);
            mult = pagemtimes(dNdxT, mult);
            intI = mult.*dV;
            lhsC = squeezeParticular(sum(intI,3),3);
%             lhsC = reshape(lhsC, [nDimf*nNodeE, nElem]);
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
