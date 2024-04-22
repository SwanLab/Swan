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

            dNdxTrial = obj.arrangedNdx(obj.trial);
            dNdxTest = obj.arrangedNdx(obj.test);
            % dNdxTest  = obj.test.evaluateCartesianDerivatives(obj.quadrature.posgp);
            % dNdxTrial = obj.trial.evaluateCartesianDerivatives(obj.quadrature.posgp);

            % dNdxTr = full(sparse(1:numel(dNdxTrial), repmat(1:size(dNdxTrial,1),1,size(dNdxTrial,2)), dNdxTrial(:)))';
            % dNdxTe = full(sparse(1:numel(dNdxTest), repmat(1:size(dNdxTest,1),1,size(dNdxTest,2)), dNdxTest(:)))';

            dNdxT = permute(dNdxTest, [2 1 3 4]);

            secPiola = permute(secPiola, [2 1 3 4]);
            mult = pagemtimes(secPiola, dNdxTrial);
            mult = pagemtimes(dNdxT, mult); % hmmm
            intI = mult.*dV;
            lhsC = squeezeParticular(sum(intI,3),3);
%             lhsC = reshape(lhsC, [nDimf*nNodeE, nElem]);
            % rhsC = zeros(nNode*nDim,nElem);
        end
       
        function B = arrangedNdx(obj, fun)
            dNdx= fun.evaluateCartesianDerivatives(obj.quadrature.posgp);
            nDimG = size(dNdx,1);
            nNode = size(dNdx,2);
            nGaus = size(dNdx,3);
            nElem = size(dNdx,4);
            B = zeros(3,24,nGaus,nElem);
            
            for iNode = 1:nNode
                add = 3*(iNode-1);
                B(1,add+1,:,:) = dNdx(1,iNode,:,:);
                B(2,add+2,:,:) = dNdx(2,iNode,:,:);
                B(3,add+3,:,:) = dNdx(3,iNode,:,:);
            end

        end

    end
    
end
