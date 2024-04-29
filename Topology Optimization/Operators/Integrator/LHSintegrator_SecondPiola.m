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
            obj.lambda = cParams.material.lambda;
            obj.mu     = cParams.material.mu;
        end
        
        function lhsC = computeElementalLHS(obj)
%             obj.mu = 10 / (2*(1 + 0.3));
%             obj.lambda = (10*0.3) / ((1 + 0.3) * (1 - 2*0.3));

            % Matrices
            uFun = obj.trial;
            quad = obj.quadrature;

            xG = quad.posgp;
            dV(1,1,:,:,:,:) = obj.mesh.computeDvolume(obj.quadrature);
            dNdxTrial(1,1,:,:,:,:) = obj.arrangedNdx(obj.trial);
            dNdxTest(1,1,:,:,:,:)  = obj.arrangedNdx(obj.test); % tesi ester 2 . 2 2 25, holzapfel te ctan

            nPoints  = obj.quadrature.ngaus;
            nElem = obj.mesh.nelem;
            nDimG = obj.mesh.ndim;
            nDimf = uFun.ndimf;
            GradU = reshape(Grad(uFun).evaluate(xG),[nDimG,nDimf,nPoints, nElem]);
%             GradU = [0.0 0.0 0.0; -3.415063509461096 -0.24999999999999956 -0.4330127018922192; 0.9150635094610968 0.43301270189221924 -0.24999999999999994];

            I33 = zeros(size(GradU));
            I33(1,1,:,:) = 1;
            I33(2,2,:,:) = 1;

            I33(3,3,:,:) = 1;
            F = I33 + GradU; % deformation gradient
            Ft = permute(F, [2 1 3 4]);
            C = pagemtimes(Ft,F);
            invC  = MatrixVectorizedInverter.computeInverse(C);
            invF = MatrixVectorizedInverter.computeInverse(F);
            invFt = MatrixVectorizedInverter.computeInverse(Ft);
            jac(1,1,:,:)  = MatrixVectorizedInverter.computeDeterminant(F);
            logJac(1,1,:,:,:,:) = log(jac);
            Aneo = obj.lambda*obj.outerProduct(invFt, invFt) + ...
                obj.mu*obj.kron_topF(I33,I33) + ...
                (obj.mu-obj.lambda*logJac).*obj.kron_botF(invFt, invF);

            % dNdxTest  = obj.test.evaluateCartesianDerivatives(obj.quadrature.posgp);
            % dNdxTrial = obj.trial.evaluateCartesianDerivatives(obj.quadrature.posgp);

            dNdxT = permute(dNdxTest, [1 2 4 3 5 6 ]);
            part = zeros([24,24,size(Aneo,5),size(Aneo,6)])
            for i = 1:size(Aneo,1)
                for j = 1:size(Aneo,2)
                    Ctan = squeeze(Aneo(i,j,:,:,:,:)); % 3x3xnGxnElem
                    pr = squeeze(sum(pagemtimes(Ctan, squeeze(dNdxTr)),1));
                    
                end
            end

            % dNdxTr = full(sparse(1:numel(dNdxTrial), repmat(1:size(dNdxTrial,1),1,size(dNdxTrial,2)), dNdxTrial(:)))';
            % dNdxTe = full(sparse(1:numel(dNdxTest), repmat(1:size(dNdxTest,1),1,size(dNdxTest,2)), dNdxTest(:)))';


            % secPiola = permute(secPiola, [2 1 3 4]);
            mult = pagemtimes(Aneo, dNdxTrial);
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

        function C = outerProduct(obj, A, B)  % version 1
            C = zeros([size(A,1), size(A,2),size(B)]);
            for i = 1:size(A,1)
                for j = 1:size(A,2)
                    for k = 1:size(B,1)
                        for l = 1:size(B,2)
                            C(i,j,k,l,:,:) = A(i,j,:,:).*B(k,l,:,:);
                        end
                    end
                end
            end
        end

        function C= kron_topF(obj, A,B)
            C = zeros([size(A,1), size(A,2), size(B)]); % to support 4th order tensors
            for i = 1:size(A,1)
                for k = 1:size(B,1)
                    C(i,:,k,:,:,:) = pagemtimes( A(i,k,:,:), B);
                end
            end
        end

        function C= kron_botF(obj,A,B)
%             C = obj.kron_topF(A,B);
%             C = pagetranspose(C);
        C = zeros([size(A,1), size(A,2), size(B)]); % to support 4th order tensors
        for i = 1:size(A,1)
            for j = 1:size(A,2)
                for k = 1:size(B,1)
                    for l = 1:size(B,2)
                        C(i,j,k,l,:,:) = A(i,l,:,:).*B(j,k,:,:);
                    end
                end
            end
        end
        end


    end
    
end
