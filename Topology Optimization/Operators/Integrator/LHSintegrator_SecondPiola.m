classdef LHSintegrator_SecondPiola < LHSintegrator
    
    % \int Piola : Fdot

    properties (Access = private)
        mu, lambda
    end

    methods (Access = public)
        
        function obj = LHSintegrator_SecondPiola(cParams)
            obj@LHSintegrator(cParams)
            obj.initMat(cParams)
        end
        
        function LHS = compute(obj)
            lhs = obj.computeElementalLHS();
            LHS = obj.assembleMatrix(lhs);
        end
        
    end
    
    methods (Access = private)
        
        function initMat(obj, cParams)
            obj.lambda = 1;
            obj.mu = 1;
        end
        
        function lhsC = computeElementalLHS(obj)
            % Matrices
            uFun = obj.trial;
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
            dNdxTest  = obj.test.evaluateCartesianDerivatives(obj.quadrature.posgp);
            dNdxTrial = obj.trial.evaluateCartesianDerivatives(obj.quadrature.posgp);

            dNdxT = permute(dNdxTest, [2 1 3 4]);
            secPiola = permute(secPiola, [2 1 3 4]);
            mult = pagemtimes(secPiola, dNdxTrial);
            mult = pagemtimes(dNdxT, mult);
            intI = mult.*dV;
            lhsC = squeezeParticular(sum(intI,3),3);
%             lhsC = reshape(lhsC, [nDimf*nNodeE, nElem]);
            % rhsC = zeros(nNode*nDim,nElem);
        end
        
    end
    
end
