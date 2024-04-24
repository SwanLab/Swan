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
            obj.lambda = cParams.material.lambda;
            obj.mu     = cParams.material.mu;
            obj.quadratureOrder = 1;
        end
        
        function rhsC = computeElementalRHS(obj, uFun, test)
            % Matrices
            quad = obj.quadrature;
%             uFun = smallU;
            xG = quad.posgp;
            nPoints  = obj.quadrature.ngaus;
            nElem = obj.mesh.nelem;
            nDimG = obj.mesh.ndim;
            nDimf = uFun.ndimf;
            GradU = reshape(Grad(uFun).evaluate(xG),[nDimG,nDimf,nPoints, nElem]);

            I33 = zeros(size(GradU));
            I33(1,1,:,:) = 1/nPoints;
            I33(2,2,:,:) = 1/nPoints;
            I33(3,3,:,:) = 1/nPoints;
            
            F = I33 + GradU; % deformation gradient
            invF = MatrixVectorizedInverter.computeInverse(F);
            invFt = permute(invF, [2 1 3 4]);
            jac(1,1,:,:)  = MatrixVectorizedInverter.computeDeterminant(F);
            
            piola = obj.mu*(F-invFt) + obj.lambda*(jac-1).*jac.*invFt;
%             piola = obj.mu*(F-invFt) + obj.lambda*(jac-1).*invFt;
            % can i subtract 1 to the jacobian directly?

%             strs = obj.revertVoigtNotation(smallStress);
%             Ft = permute(F, [2 1 3 4]);
%             streMeu = pagemtimes(1./jac, pagemtimes(Ft,piola))
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

    methods (Static, Access = private)

        function f_new = revertVoigtNotation(f)
            ndim = size(f,1);
            ngaus = size(f,2);
            nel  = size(f,3);
            b = 1;
            switch ndim
                case 1
                    f_new = f;
                case 3
                    f_new        = zeros(2,2,ngaus,nel);
                    f_new(1,1,:,:) = f(1,:,:,:);
                    f_new(2,2,:,:) = f(2,:,:,:);
                    f_new(1,2,:,:) = f(3,:,:,:)/b;
                    f_new(2,1,:,:) = f_new(1,2,:,:);
                case 6
                    f_new        = zeros(3,3,ngaus,nel);
                    f_new(1,1,:,:) = f(1,:,:,:);
                    f_new(2,2,:,:) = f(2,:,:,:);
                    f_new(3,3,:,:) = f(3,:,:,:);
                    f_new(1,2,:,:) = f(6,:,:,:)/b;
                    f_new(2,1,:,:) = f_new(1,2,:,:);
                    f_new(1,3,:,:) = f(5,:,:,:)/b;
                    f_new(3,1,:,:) = f_new(1,3,:,:);
                    f_new(2,3,:,:) = f(4,:,:,:)/b;
                    f_new(3,2,:,:) = f_new(2,3,:,:);
            end
        end
        
    end
    
end
