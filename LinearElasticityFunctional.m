classdef LinearElasticityFunctional < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        mesh
        lambda
        mu
        material
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = LinearElasticityFunctional (cParams)
            obj.init(cParams)
%             obj.mu = 10 / (2*(1 + 0.3));
%             obj.lambda = (10*0.3) / ((1 + 0.3) * (1 - 2*0.3));
        end

        function val = compute(obj, uFun)
            quad = Quadrature.create(obj.mesh, 2);
            xG = quad.posgp;
            nDimf = uFun.ndimf;
            
            [F,~] = obj.computeDeformationGradient(uFun, xG);
            Ft = permute(F, [2 1 3 4]);
            
            C = pagemtimes(Ft,F);
            trC = obj.computeTrace(C);

            jac(1,1,:,:)  = MatrixVectorizedInverter.computeDeterminant(F);
            dV(1,1,:,:) = obj.mesh.computeDvolume(quad);

            % val = obj.mu/2*(trC - nDimf) - obj.mu*log(jac) + obj.lambda/2*(jac-1).^2;
            val = obj.mu/2*(trC - nDimf) - obj.mu*log(jac) + obj.lambda/2*(log(jac)).^2; % stanford
            val = pagemtimes(val,dV);
            val = sum(val, 'all');
        end
        

        function Ju = computeGradient(obj, uFun)
%             quad = Quadrature.create(obj.mesh, 2);
%             xV = quad.posgp;
%             C = obj.material.obtainTensor();
            sigma = DDP(obj.material,Voigt(SymGrad(uFun)));
            test = LagrangianFunction.create(obj.mesh, uFun.ndimf, uFun.order);

            s.mesh = obj.mesh;
            s.quadratureOrder = 1;
            s.type = 'ShapeSymmetricDerivative';
            RHS = RHSintegrator.create(s);
            Ju = RHS.compute(sigma,test);
        end

        function Huu = computeHessian(obj, uFun) 
            s.type     = 'ElasticStiffnessMatrix';
            s.mesh     = obj.mesh;
            s.material = obj.material;
            s.quadratureOrder = 2;
            s.test     = LagrangianFunction.create(obj.mesh,uFun.ndimf, 'P1');
            s.trial    = uFun;
            LHS = LHSintegrator.create(s);
            Huu = LHS.compute();
        end


    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.lambda = cParams.lambda;
            obj.mu     = cParams.mu;
            obj.mesh   = cParams.mesh;
            obj.material = cParams.material;
        end

        function piola = computeFirstPiola(obj,uFun,xG)
            [F,~] = obj.computeDeformationGradient(uFun, xG);
            Ft = permute(F, [2 1 3 4]);

            nGaus = size(xG,2);
            nElem = size(Ft,4);

            invFt = MatrixVectorizedInverter.computeInverse(Ft);
            jac(1,1,:,:)  = MatrixVectorizedInverter.computeDeterminant(F);
            jac = reshape(jac, [1 1 nGaus nElem]);
            piola = obj.mu*(F-invFt) + obj.lambda*log(jac).*invFt;
%             piola = obj.mu*(F-invFt) + obj.lambda/2*(jac.^2-1).*invFt;
        end

        function Aneo = computeTangentConstitutive(obj,uFun,xG)
            [F,I33] = obj.computeDeformationGradient(uFun, xG);
            Ft = permute(F, [2 1 3 4]);
            
            nGaus = size(xG,2);
            nElem = size(Ft,4);

            invF = MatrixVectorizedInverter.computeInverse(F);
            invFt = MatrixVectorizedInverter.computeInverse(Ft);

            jac(1,1,:,:)  = MatrixVectorizedInverter.computeDeterminant(F);
            jac = reshape(jac, [1 1 1 1 nGaus nElem]);
            logJac = log(jac);

            Aneo = obj.lambda*obj.outerProduct(invFt, invFt) + ...
                obj.mu*obj.kron_topF(I33,I33) + ...
                (obj.mu-obj.lambda*logJac).*obj.kron_botF(invFt, invF);
%             quad = Quadrature.create(obj.mesh,3);
%             xG = quad.posgp;

%             [F,~] = obj.computeDeformationGradient(uFun, xG);
%             Ft = permute(F, [2 1 3 4]);
%             nF    = size(F,1);
%             nGaus = size(F,3);
%             nElem = size(F,4);
%             A = zeros(nF,nF,nF,nF,nGaus,nElem);
%                 for i = 1:nF
%                 for j = 1:nF
%                 for k = 1:nF
%                 for l = 1:nF
%                     first_term(1,1,1,1,:,:) = squeeze(obj.lambda*invF(j,i,:,:).*invF(l,k,:,:));
%                     secnd_term(1,1,1,1,:,:) = squeeze(obj.mu*I33(i,k,:,:).*I33(j,l,:,:));
%                     trd(1,1,:,:,:,:) = invF(l,i,:,:).*invF(j,k,:,:);
%                     third_term(1,1,1,1,:,:) = squeeze((obj.mu - obj.lambda*log(jac)).*trd);
%                     A(i,j,k,l,:,:) = first_term + secnd_term + third_term;
% %                     A(i,j,k,l,:,:) = first_term;
%                 end
%                 end
%                 end
%                 end
%             norm(Aneo(:)-A(:))
%             Aneo = A;
        end

        function [F,I33] = computeDeformationGradient(obj, uFun, xG)
            nPoints  = size(xG,2);
            nElem = obj.mesh.nelem;
            nDimG = obj.mesh.ndim;
            nDimf = uFun.ndimf;

            GradU = reshape(Grad(uFun).evaluate(xG),[nDimG,nDimf,nPoints, nElem]);
            GradU = permute(GradU, [2 1 3 4]);
%             GradU = [0.0 0.0 0.0; -3.415063509461096 -0.24999999999999956 -0.4330127018922192; 0.9150635094610968 0.43301270189221924 -0.24999999999999994];

            I33 = obj.createIdentityMatrix(size(GradU));

            F = I33 + GradU;
        end

        function trC = computeTrace(obj,C)
            trC = zeros(1,1,size(C,3),size(C,4));
            for i = 1:size(C,1)
                trC = trC + C(i,i,:,:,:);
            end
        end

        function I = createIdentityMatrix(obj,sz)
            I = zeros(sz);
            for i = 1:sz(1)
                I(i,i,:,:) = 1;
            end
        end

        % Operators
        function C = double_dot(obj,A, B)
            assert(~isvector(A) && ~isvector(B))
            idx = max(0, ndims(A) - 1);
            B_t = permute(B, circshift(1:ndims(A) + ndims(B), [0, idx - 1]));
            C = squeeze(sum(squeeze(sum(bsxfun(@times, A, B_t), idx)), idx));
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

%         function C= kron_topF(obj, A,B)
%             C = zeros([size(A,1), size(A,2), size(B)]); % to support 4th order tensors
%             for i = 1:size(A,1)
%                 for k = 1:size(B,1)
%                     C(i,:,k,:,:,:) = pagemtimes( A(i,k,:,:), B);
%                 end
%             end
%         end

        function C= kron_topF(obj,A,B)
    %             C = obj.kron_topF(A,B);
    %             C = pagetranspose(C);
            C = zeros([size(A,1), size(A,2), size(B)]); % to support 4th order tensors
            for i = 1:size(A,1)
                for k = 1:size(A,2)
                    for j = 1:size(B,1)
                        for l = 1:size(B,2)
                            C(i,j,k,l,:,:) = A(i,k,:,:).*B(j,l,:,:);
                        end
                    end
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